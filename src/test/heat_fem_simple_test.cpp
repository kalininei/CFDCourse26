#include "cfd/debug/printer.hpp"
#include "cfd/debug/saver.hpp"
#include "cfd/fem/fem_assembler.hpp"
#include "cfd/fem/fem_cross_assembler.hpp"
#include "cfd/numeric_integration/triangle_quadrature.hpp"
#include "cfd/fem/elem2d/triangle_linear.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/umfpack_solver.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd/mat/matrix_iter.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd26_test.hpp"
#include "utils/filesystem.hpp"
#include "cfd/mat/algebraic_bc.hpp"
#include "utils/vecmat.hpp"
#include <list>

using namespace cfd;

namespace {

constexpr double ALPHA_U = 0.8;
constexpr double ALPHA_P = 0.3;

struct HeatFemSimpleWorker {
    HeatFemSimpleWorker(const IGrid& grid, double Re, double Pe, double delta_t);
    void initialize_saver(std::string stem, double timestep);
    void init_timestep();
    double init_step();
    void step();
    void save_current_fields(double time);
private:
    const IGrid& grid_;
    const double Re_;
    const double Pe_;
    const double delta_t_;
    FemAssembler fem_a_;
    FemAssembler fem_b_;
    FemCrossAssembler fem_cross_;

    std::vector<size_t> heater_nodes_, cooler_nodes_;
    CsrMatrix A_, B_, Bt_, Z_, H_, Schur_, HBt_;
    std::vector<double> rhs_;
    std::vector<double> uvp_;
    std::vector<double> uvp_old_;
    std::vector<double> temperature_;
    std::vector<double> temperature_old_;

    double* u_;
    double* v_;
    double* p_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_;

    void assemble_slae();
    void solve_temperature_problem();
    CsrMatrix assemble_A() const;
    CsrMatrix assemble_B() const;
    std::vector<double> assemble_rhs() const;
    std::vector<double> use_preconditioner(const std::vector<double>& r) const;
    std::vector<double> compute_residual() const;

    static FemAssembler build_fem_a(const IGrid& grid);
    static FemAssembler build_fem_b(const IGrid& grid);
};

HeatFemSimpleWorker::HeatFemSimpleWorker(const IGrid& grid, double Re, double Pe, double delta_t)
    : grid_(grid),
      Re_(Re),
      Pe_(Pe),
      delta_t_(delta_t),
      fem_a_(build_fem_a(grid_)),
      fem_b_(build_fem_b(grid_)),
      fem_cross_(fem_b_, fem_a_){

    // boundary nodes
    double ymax = grid_.box().second.y;
    for (size_t inode: grid_.boundary_points()){
        Point p = grid_.point(inode);
        double r = vector_abs(p);
        if (p.y > ymax - 0.0001){
            cooler_nodes_.push_back(inode);
        } else if (r < 1.001){
            heater_nodes_.push_back(inode);
        }
    }

    // init vectors
    uvp_ = std::vector<double>(fem_a_.n_bases() + fem_b_.n_bases(), 0);
    rhs_ = std::vector<double>(uvp_.size());
    u_ = uvp_.data();
    v_ = u_ + fem_a_.n_bases()/2;
    p_ = u_ + fem_a_.n_bases();
    temperature_.resize(fem_b_.n_bases(), 0);

    // ============= Matrices
    // Z
    LodMatrix zero(grid_.n_points());
    for (size_t i=0; i<grid_.n_points(); ++i){
        zero.set_value(i, i, 0.0);
    }
    Z_ = zero.to_csr();
    // B
    B_ = assemble_B();
    // dirichlet boundary conditions
    double* rhs_b_ = rhs_.data() + 2*(grid_.n_points() + grid_.n_cells());
    for (size_t i: grid_.boundary_points()){
        size_t u_idx = i;
        size_t v_idx = grid_.n_points() + grid_.n_cells() + i;
        algebraic_bc_exclude_column(u_idx, 0.0, B_, rhs_b_);  // since u = 0 we dont change rhs_b here
        algebraic_bc_exclude_column(v_idx, 0.0, B_, rhs_b_);
    }
    // p[0] = 0
    B_.clear_row(0);
    Z_.set_value(0, 0, 1.0);
    // Bt
    Bt_ = mat_transpose(B_);
}

FemAssembler HeatFemSimpleWorker::build_fem_a(const IGrid& grid) {
    constexpr int GAUSS_POWER = 3;
    const size_t iu0 = 0;
    const size_t iv0 = grid.n_points() + grid.n_cells();
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        std::vector<size_t> ipoints = grid.tab_cell_point(icell);
        FemElement el;
        el.geometry = std::make_shared<TriangleLinearGeometry>(grid.point(ipoints[0]), grid.point(ipoints[1]), grid.point(ipoints[2]));
        std::vector<std::shared_ptr<IElementBasis>> basis_list{
            std::make_shared<TriangleLinearBubbleBasis>(),  // u
            std::make_shared<TriangleLinearBubbleBasis>()}; // v
        el.basis = std::make_shared<CompaundBasis>(basis_list);
        el.quadrature = quadrature_triangle_gauss<GAUSS_POWER>();
        elements.push_back(el);

        std::vector<size_t> elem_bases = {
            iu0 + ipoints[0], iu0 + ipoints[1], iu0 + ipoints[2], iu0 + grid.n_points() + icell,  // u
            iv0 + ipoints[0], iv0 + ipoints[1], iv0 + ipoints[2], iv0 + grid.n_points() + icell   // v
        };
        tab_elem_basis.push_back(elem_bases);
    }
    size_t n_bases = 2 * grid.n_points() + 2 * grid.n_cells();
    return FemAssembler(n_bases, elements, tab_elem_basis);
}

FemAssembler HeatFemSimpleWorker::build_fem_b(const IGrid& grid) {
    constexpr int GAUSS_POWER = 3;
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        std::vector<size_t> ipoints = grid.tab_cell_point(icell);
        FemElement el;
        el.geometry = std::make_shared<TriangleLinearGeometry>(grid.point(ipoints[0]), grid.point(ipoints[1]), grid.point(ipoints[2]));
        el.basis = std::make_shared<TriangleLinearBasis>();
        el.quadrature = quadrature_triangle_gauss<GAUSS_POWER>();
        elements.push_back(el);

        std::vector<size_t> elem_bases = ipoints;
        tab_elem_basis.push_back(elem_bases);
    }
    size_t n_bases = grid.n_points();
    return FemAssembler(n_bases, elements, tab_elem_basis);
}

void HeatFemSimpleWorker::initialize_saver(std::string stem, double timestep) {
    writer_.reset(new VtkUtils::TimeSeriesWriter(stem));
    writer_->set_time_step(timestep);
};

CsrMatrix HeatFemSimpleWorker::assemble_A() const{
    CsrMatrix ret;
    ret.set_stencil(fem_a_.stencil());

    for (size_t ielem = 0; ielem < fem_a_.n_elements(); ++ielem) {
        const FemElement& el = fem_a_.element(ielem);
        std::vector<double> local_uv_old = fem_a_.local_vector(ielem, uvp_);

        auto fun = [&](Point p) -> std::vector<double> {
            const size_t n = el.basis->size();
            std::vector<double> loc(n * n, 0.0);
            auto val = FemElementValue(&el, p);
            // convection velocity
            const Vector a = {val.subrange_interpolate(0, 4, local_uv_old),
                              val.subrange_interpolate(4, 8, local_uv_old)};
            // u
            for (size_t ibas = 0; ibas < 4; ++ibas) {
                Vector g1 = val.grad_phi(ibas);
                // u/dt (lumped)
                loc[ibas * n + ibas] += val.phi(ibas) * val.modj() / delta_t_;
                for (size_t jbas = 0; jbas < 4; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // -1/Re * laplace(u)
                    Vector g2 = val.grad_phi(jbas);
                    loc[k1] += (1.0/Re_) * dot_product(g1, g2) * val.modj();
                    // a * grad(u)
                    loc[k1] += dot_product(a, g2) * val.phi(ibas) * val.modj();
                }
            }
            // v
            for (size_t ibas = 4; ibas < 8; ++ibas) {
                Vector g1 = val.grad_phi(ibas);
                // v/dt (lumped)
                loc[ibas * n + ibas] += val.phi(ibas) * val.modj() / delta_t_;
                for (size_t jbas = 4; jbas < 8; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // -1/Re * laplace(v)
                    Vector g2 = val.grad_phi(jbas);
                    loc[k1] += (1.0/Re_) * dot_product(g1, g2) * val.modj();
                    // a * grad(v)
                    loc[k1] += dot_product(a, g2) * val.phi(ibas) * val.modj();
                }
            }
            return loc;
        };
        auto local_mat = el.quadrature->integrate(fun);
        fem_a_.add_to_global_matrix(ielem, local_mat, ret.vals());
    }

    return ret;
}

CsrMatrix HeatFemSimpleWorker::assemble_B() const {
    CsrMatrix ret;
    ret.set_stencil(fem_cross_.stencil());

    for (size_t ielem = 0; ielem < fem_b_.n_elements(); ++ielem) {
        const FemElement& el_a = fem_a_.element(ielem);
        const FemElement& el_b = fem_b_.element(ielem);

        auto fun = [&](Point p) -> std::vector<double> {
            const size_t na = el_a.basis->size();
            const size_t nb = el_b.basis->size();
            std::vector<double> loc(na * nb, 0.0);
            auto val_a = FemElementValue(&el_a, p);
            auto val_b = FemElementValue(&el_b, p);

            for (size_t ibas = 0; ibas < el_b.basis->size(); ++ibas) {
                // du/dx
                for (size_t jbas=0; jbas < na/2; ++jbas){
                    const size_t k1 = ibas * na + jbas;
                    Vector g2 = val_a.grad_phi(jbas);
                    loc[k1] += g2.x * val_b.phi(ibas) * val_b.modj();
                }
                // dv/dx
                for (size_t jbas=na/2; jbas < na; ++jbas){
                    const size_t k1 = ibas * na + jbas;
                    Vector g2 = val_a.grad_phi(jbas);
                    loc[k1] += g2.y * val_b.phi(ibas) * val_b.modj();
                }
            }
            return loc;
        };
        auto local_mat = el_b.quadrature->integrate(fun);
        fem_cross_.add_to_global_matrix(1.0, ielem, local_mat, ret.vals());
    }

    return ret;
}

std::vector<double> HeatFemSimpleWorker::assemble_rhs() const{
    std::vector<double> ret(uvp_.size(), 0);

    for (size_t ielem = 0; ielem < fem_a_.n_elements(); ++ielem) {
        const FemElement& el_a = fem_a_.element(ielem);
        const FemElement& el_b = fem_b_.element(ielem);
        std::vector<double> local_uv_old = fem_a_.local_vector(ielem, uvp_);
        std::vector<double> local_temp = fem_b_.local_vector(ielem, temperature_);

        auto fun = [&](Point p) -> std::vector<double> {
            const size_t n = el_a.basis->size();
            std::vector<double> loc(n, 0.0);
            auto val_a = FemElementValue(&el_a, p);
            auto val_b = FemElementValue(&el_b, p);
            // temperature
            const double temp = val_b.interpolate(local_temp);
            // velocity old
            const double uold = val_a.subrange_interpolate(0, 4, local_uv_old);
            const double vold = val_a.subrange_interpolate(4, 8, local_uv_old);
            // u
            for (size_t ibas = 0; ibas < 4; ++ibas) {
                // Uold/tau
                loc[ibas] += (uold/delta_t_) * val_a.phi(ibas) * val_a.modj();
            }
            // v
            for (size_t ibas = 4; ibas < 8; ++ibas) {
                // Vold/tau + T
                loc[ibas] += (vold/delta_t_ + temp) * val_a.phi(ibas) * val_a.modj();
            }
            return loc;
        };
        auto local_vec = el_a.quadrature->integrate(fun);
        fem_a_.add_to_global_vector(ielem, local_vec, ret);
    }

    return ret;
}

std::vector<double> HeatFemSimpleWorker::use_preconditioner(const std::vector<double>& r) const{
    CsrMatrix B01 = mat_multiply(A_, HBt_);
    CsrMatrix precond = assemble_block_matrix({{&A_, &B01}, {&B_, &Z_}});
    std::vector<double> x(r.size(), 0);
    UmfpackSolver::solve_slae(precond, r, x);

    return x;
}

void HeatFemSimpleWorker::init_timestep() {
    temperature_old_ = temperature_;
    uvp_old_ = uvp_;
    solve_temperature_problem();
}

void HeatFemSimpleWorker::solve_temperature_problem(){
    // assembly
    CsrMatrix mat;
    mat.set_stencil(fem_b_.stencil());
    std::vector<double> rhs(fem_b_.n_bases(), 0);

    for (size_t ielem = 0; ielem < fem_b_.n_elements(); ++ielem) {
        const FemElement& el_a = fem_a_.element(ielem);
        const FemElement& el_b = fem_b_.element(ielem);
        std::vector<double> local_uv_old = fem_a_.local_vector(ielem, uvp_old_);
        std::vector<double> local_temp_old = fem_b_.local_vector(ielem, temperature_old_);

        // LHS
        auto fun_lhs = [&](Point p) -> std::vector<double> {
            const size_t n = el_b.basis->size();
            std::vector<double> loc(n*n, 0.0);
            auto val_a = FemElementValue(&el_a, p);
            auto val_b = FemElementValue(&el_b, p);
            // convection velocity
            const Vector a = {val_a.subrange_interpolate(0, 4, local_uv_old),
                              val_a.subrange_interpolate(4, 8, local_uv_old)};
            for (size_t ibas = 0; ibas < n; ++ibas) {
                Vector g1 = val_b.grad_phi(ibas);
                // t/dt (lumped)
                loc[ibas*n + ibas] += val_b.phi(ibas) * val_b.modj() / delta_t_;
                for (size_t jbas = 0; jbas < n; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // -1/Pe * laplace(t)
                    Vector g2 = val_b.grad_phi(jbas);
                    loc[k1] += (1.0/Pe_) * dot_product(g1, g2) * val_b.modj();
                    // a * grad(u)
                    loc[k1] += dot_product(a, g2) * val_b.phi(ibas) * val_b.modj();
                }
            }
            return loc;
        };
        auto local_mat = el_b.quadrature->integrate(fun_lhs);
        fem_b_.add_to_global_matrix(ielem, local_mat, mat.vals());

        // RHS
        auto fun_rhs = [&](Point p) -> std::vector<double> {
            const size_t n = el_b.basis->size();
            std::vector<double> loc(n, 0.0);
            auto val_b = FemElementValue(&el_b, p);
            double temp = val_b.interpolate(local_temp_old);
            for (size_t ibas = 0; ibas < n; ++ibas) {
                // t/dt
                loc[ibas] += temp * val_b.phi(ibas) * val_b.modj() / delta_t_;
            }
            return loc;
        };
        auto local_vec = el_b.quadrature->integrate(fun_rhs);
        fem_b_.add_to_global_vector(ielem, local_vec, rhs);
    }

    // boundary conditions
    for (auto inode: heater_nodes_){
        mat.set_unit_row(inode);
        rhs[inode] = 1;
    }
    for (auto inode: cooler_nodes_){
        mat.set_unit_row(inode);
        rhs[inode] = 0;
    }
    // solve slae
    AmgcMatrixSolver::solve_slae(mat, rhs, temperature_);
}

void HeatFemSimpleWorker::assemble_slae(){
    // A and rhs
    A_ = assemble_A();
    rhs_ = assemble_rhs();

    // relax A
    for (size_t irow = 0; irow < fem_a_.n_bases(); ++irow){
        double d = A_.value(irow, irow);
        A_.set_value(irow, irow, d/ALPHA_U);
        rhs_[irow] += d / ALPHA_U * (1 - ALPHA_U) * uvp_[irow];
    }

    // dirichlet boundary conditions
    for (size_t i: grid_.boundary_points()){
        size_t u_idx = i;
        size_t v_idx = grid_.n_points() + grid_.n_cells() + i;
        // we need symmetry to keep B and Bt
        symmetric_algebraic_bc_dirichlet({{u_idx, 0}, {v_idx, 0}}, A_, rhs_);
    }
    
    // H = diag(A)^{-1}   <= SIMPLE preconditioner
    LodMatrix Hlod(fem_a_.n_bases());
    for (size_t i=0; i<Hlod.n_rows(); ++i){
        Hlod.set_value(i, i, 1.0 / A_.value(i, i));
    }
    H_ = Hlod.to_csr();

    // H*Bt
    HBt_ = mat_multiply(H_, Bt_);

    // Schur = B * H * Bt
    Schur_ = mat_multiply(B_, HBt_);
    Schur_.set_value(0, 0, Schur_.value(0, 0) - 1.0);  // <- p[0] = 0
}

double HeatFemSimpleWorker::init_step() {
    assemble_slae();

    // return max(|residual|)
    double rmax = 0;
    for (double r: compute_residual()){
        rmax = std::max(rmax, std::abs(r));
    }
    return rmax;
}

std::vector<double> HeatFemSimpleWorker::compute_residual() const{
    CsrMatrix mat = assemble_block_matrix({{&A_, &Bt_}, {&B_, &Z_}});
    // compute residual: rhs_ - M*uvp
    std::vector<double> residual = rhs_;
    for (auto [i, j, aij]: matrix_iter::ijv(mat)){
        residual[i] -= aij * uvp_[j];
    }
    return residual;
}

void HeatFemSimpleWorker::step() {
    std::vector<double> residual = compute_residual();

    // compute correction
    std::vector<double> delta_uvp = use_preconditioner(residual);

    // use correction
    for (size_t i=0; i<uvp_.size(); ++i){
        if (i < fem_a_.n_bases()){
            // u correction
            uvp_[i] += delta_uvp[i];
        } else{
            // p correction
            uvp_[i] += ALPHA_P * delta_uvp[i];
        }
    }
}

void HeatFemSimpleWorker::save_current_fields(double tm) {
    if (writer_) {
        std::string filepath = writer_->add(tm);
        grid_.save_vtk(filepath);
        const size_t np = grid_.n_points();
        VtkUtils::add_point_data(std::vector<double>(p_, p_ + np), "pressure", filepath, grid_.n_points());
        VtkUtils::add_point_data(temperature_, "temperature", filepath, grid_.n_points());
        VtkUtils::add_point_vector(std::vector<double>(u_, u_ + np), std::vector<double>(v_, v_ + np), "velocity",
                                  filepath, grid_.n_points());
    }
}

} // namespace

TEST_CASE("Heat FEM-SIMPLE", "[heat-fem-simple]") {
    std::cout << std::endl << "--- cfd26_test [heat-fem-simple] --- " << std::endl;

    // problem parameters
    double Re = 100;
    double Pe = 100;
    size_t max_it = 10'000;
    double eps = 1e-3;
    double dt = 1e-1;
    double tend = 0.3;

    // worker initialization
    const std::string gridfn = test_directory_file("cylgrid_vertical_3k.vtk");
    auto grid = UnstructuredGrid2D::vtk_read(gridfn);
    HeatFemSimpleWorker worker(grid, Re, Pe, dt);
    worker.initialize_saver("heat-fem-simple", 0.0);
    worker.save_current_fields(0.0);

    // time loop
    double tm = 0;
    double nrm = 0;
    while (tm < tend + 1e-6){
        std::cout << "----- time = " << tm << std::endl;
        worker.init_timestep();
        // iterations loop
        size_t it = 0;
        for (it = 1; it < max_it; ++it) {
            nrm = worker.init_step();
            std::cout << it-1 << " " << nrm << std::endl;
            if (it > 1 && nrm < eps){
                break;
            }
            worker.step();
        }
        worker.save_current_fields(tm);
        tm += dt;
    }
    CHECK(nrm == Approx(0.0009008417));
}
