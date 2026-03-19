#include "cfd/debug/printer.hpp"
#include "cfd/debug/saver.hpp"
#include "cfd/fem/elem2d/triangle_linear.hpp"
#include "cfd/fem/fem_assembler.hpp"
#include "cfd/fem/fem_cross_assembler.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/algebraic_bc.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/matrix_iter.hpp"
#include "cfd/mat/umfpack_solver.hpp"
#include "cfd/numeric_integration/triangle_quadrature.hpp"
#include "cfd26_test.hpp"
#include "utils/filesystem.hpp"
#include "utils/vecmat.hpp"
#include <list>

using namespace cfd;

namespace {

struct CavityFemSimpleWorker {
    CavityFemSimpleWorker(const IGrid& grid, double Re);
    void initialize_saver(std::string stem);
    double init_step();
    void step();
    void save_current_fields(size_t iter);

private:
    const IGrid& grid_;
    const double Re_;
    FemAssembler fem_a_;
    FemAssembler fem_b_;
    FemCrossAssembler fem_cross_;

    CsrMatrix A_, B_, Bt_, H_, Schur_, HBt_;
    CsrMatrix mat_;
    std::vector<double> rhs_;
    std::vector<double> uvp_;

    double* u_;
    double* v_;
    double* p_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_;

    CsrMatrix assemble_A() const;
    CsrMatrix assemble_B() const;
    std::vector<double> use_preconditioner(const std::vector<double>& r) const;

    static FemAssembler build_fem_a(const IGrid& grid);
    static FemAssembler build_fem_b(const IGrid& grid);
};

CavityFemSimpleWorker::CavityFemSimpleWorker(const IGrid& grid, double Re)
    : grid_(grid),
      Re_(Re),
      fem_a_(build_fem_a(grid_)),
      fem_b_(build_fem_b(grid_)),
      fem_cross_(fem_b_, fem_a_) {

    uvp_ = std::vector<double>(fem_a_.n_bases() + fem_b_.n_bases(), 0);
    rhs_ = std::vector<double>(uvp_.size());
    u_ = uvp_.data();
    v_ = u_ + fem_a_.n_bases() / 2;
    p_ = u_ + fem_a_.n_bases();
}

FemAssembler CavityFemSimpleWorker::build_fem_a(const IGrid& grid) {
    constexpr int GAUSS_POWER = 3;
    const size_t iu0 = 0;
    const size_t iv0 = grid.n_points() + grid.n_cells();
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        std::vector<size_t> ipoints = grid.tab_cell_point(icell);
        FemElement el;
        el.geometry = std::make_shared<TriangleLinearGeometry>(grid.point(ipoints[0]), grid.point(ipoints[1]),
                                                               grid.point(ipoints[2]));
        std::vector<std::shared_ptr<IElementBasis>> basis_list{std::make_shared<TriangleLinearBubbleBasis>(),  // u
                                                               std::make_shared<TriangleLinearBubbleBasis>()}; // v
        el.basis = std::make_shared<CompaundBasis>(basis_list);
        el.quadrature = quadrature_triangle_gauss<GAUSS_POWER>();
        elements.push_back(el);

        std::vector<size_t> elem_bases = {
            iu0 + ipoints[0], iu0 + ipoints[1], iu0 + ipoints[2], iu0 + grid.n_points() + icell, // u
            iv0 + ipoints[0], iv0 + ipoints[1], iv0 + ipoints[2], iv0 + grid.n_points() + icell  // v
        };
        tab_elem_basis.push_back(elem_bases);
    }
    size_t n_bases = 2 * grid.n_points() + 2 * grid.n_cells();
    return FemAssembler(n_bases, elements, tab_elem_basis);
}

FemAssembler CavityFemSimpleWorker::build_fem_b(const IGrid& grid) {
    constexpr int GAUSS_POWER = 3;
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        std::vector<size_t> ipoints = grid.tab_cell_point(icell);
        FemElement el;
        el.geometry = std::make_shared<TriangleLinearGeometry>(grid.point(ipoints[0]), grid.point(ipoints[1]),
                                                               grid.point(ipoints[2]));
        el.basis = std::make_shared<TriangleLinearBasis>();
        el.quadrature = quadrature_triangle_gauss<GAUSS_POWER>();
        elements.push_back(el);

        std::vector<size_t> elem_bases = ipoints;
        tab_elem_basis.push_back(elem_bases);
    }
    size_t n_bases = grid.n_points();
    return FemAssembler(n_bases, elements, tab_elem_basis);
}

void CavityFemSimpleWorker::initialize_saver(std::string stem) {
    writer_.reset(new VtkUtils::TimeSeriesWriter(stem));
};

CsrMatrix CavityFemSimpleWorker::assemble_A() const {
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
                for (size_t jbas = 0; jbas < 4; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // -1/Re * laplace(u)
                    Vector g2 = val.grad_phi(jbas);
                    loc[k1] += (1.0 / Re_) * dot_product(g1, g2) * val.modj();
                    // a * grad(u)
                    loc[k1] += dot_product(a, g2) * val.phi(ibas) * val.modj();
                }
            }
            // v
            for (size_t ibas = 4; ibas < 8; ++ibas) {
                Vector g1 = val.grad_phi(ibas);
                for (size_t jbas = 4; jbas < 8; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // -1/Re * laplace(v)
                    Vector g2 = val.grad_phi(jbas);
                    loc[k1] += (1.0 / Re_) * dot_product(g1, g2) * val.modj();
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

CsrMatrix CavityFemSimpleWorker::assemble_B() const {
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
                for (size_t jbas = 0; jbas < na / 2; ++jbas) {
                    const size_t k1 = ibas * na + jbas;
                    Vector g2 = val_a.grad_phi(jbas);
                    loc[k1] += g2.x * val_b.phi(ibas) * val_b.modj();
                }
                // dv/dx
                for (size_t jbas = na / 2; jbas < na; ++jbas) {
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

std::vector<double> CavityFemSimpleWorker::use_preconditioner(const std::vector<double>& r) const {
    LodMatrix zero(grid_.n_points());
    for (size_t i = 0; i < grid_.n_points(); ++i) {
        zero.set_value(i, i, (i == 0) ? 1.0 : 0.0);
    }
    CsrMatrix zero1 = zero.to_csr();

    CsrMatrix B01 = mat_multiply(A_, HBt_);
    CsrMatrix precond = assemble_block_matrix({{&A_, &B01}, {&B_, &zero1}});
    std::vector<double> x(r.size(), 0);
    UmfpackSolver::solve_slae(precond, r, x);

    return x;
}

double CavityFemSimpleWorker::init_step() {
    LodMatrix zero(grid_.n_points());
    for (size_t i = 0; i < grid_.n_points(); ++i) {
        zero.set_value(i, i, 0.0);
    }
    CsrMatrix zero1 = zero.to_csr();
    rhs_ = std::vector<double>(uvp_.size(), 0.0);

    // A and B
    A_ = assemble_A();
    B_ = assemble_B();

    // relax A
    double alpha_u = 0.8;
    for (size_t irow = 0; irow < fem_a_.n_bases(); ++irow) {
        double d = A_.value(irow, irow);
        A_.set_value(irow, irow, d / alpha_u);
        rhs_[irow] += d / alpha_u * (1 - alpha_u) * uvp_[irow];
    }

    // dirichlet boundary conditions
    double* rhs_b_ = rhs_.data() + 2 * (grid_.n_points() + grid_.n_cells());
    for (size_t i: grid_.boundary_points()) {
        size_t u_idx = i;
        size_t v_idx = grid_.n_points() + grid_.n_cells() + i;
        double valu = (grid_.point(i).y > 0.9999) ? 1.0 : 0.0;
        symmetric_algebraic_bc_dirichlet({{u_idx, valu}, {v_idx, 0}}, A_, rhs_);
        algebraic_bc_exclude_column(u_idx, valu, B_, rhs_b_);
        algebraic_bc_exclude_column(v_idx, 0.0, B_, rhs_b_);
    }
    B_.clear_row(0);
    zero1.set_value(0, 0, 1.0);

    // Bt
    Bt_ = mat_transpose(B_);

    // H = diag(A)^{-1}   <= SIMPLE preconditioner
    LodMatrix Hlod(fem_a_.n_bases());
    for (size_t i = 0; i < Hlod.n_rows(); ++i) {
        Hlod.set_value(i, i, 1.0 / A_.value(i, i));
    }
    H_ = Hlod.to_csr();

    // H*Bt
    HBt_ = mat_multiply(H_, Bt_);

    // Schur = B * H * Bt
    Schur_ = mat_multiply(B_, HBt_);
    Schur_.set_value(0, 0, Schur_.value(0, 0) - 1.0); // <- p[0] = 0

    mat_ = assemble_block_matrix({{&A_, &Bt_}, {&B_, &zero1}});
    return compute_residual(mat_, rhs_, uvp_);
}

void CavityFemSimpleWorker::step() {
    // compute residual
    std::vector<double> residual = rhs_;
    for (auto [i, j, aij]: matrix_iter::ijv(mat_)) {
        residual[i] -= aij * uvp_[j];
    }

    // compute correction
    std::vector<double> delta_uvp = use_preconditioner(residual);

    // use correction
    for (size_t i = 0; i < uvp_.size(); ++i) {
        if (i < fem_a_.n_bases()) {
            uvp_[i] += delta_uvp[i];
        } else {
            uvp_[i] += 0.3 * delta_uvp[i];
        }
    }
}

void CavityFemSimpleWorker::save_current_fields(size_t iter) {
    if (writer_) {
        std::string filepath = writer_->add_iter(iter);
        grid_.save_vtk(filepath);
        const size_t np = grid_.n_points();
        VtkUtils::add_point_data(std::vector<double>(p_, p_ + np), "pressure", filepath, grid_.n_points());
        VtkUtils::add_point_vector(std::vector<double>(u_, u_ + np), std::vector<double>(v_, v_ + np), "velocity",
                                   filepath, grid_.n_points());
    }
}

} // namespace

TEST_CASE("Cavity FEM-SIMPLE", "[cavity-fem-simple]") {
    std::cout << std::endl << "--- cfd26_test [cavity-fem-simple] --- " << std::endl;

    // problem parameters
    double Re = 100;
    size_t max_it = 10'000;
    double eps = 1e-4;

    // worker initialization
    const std::string gridfn = test_directory_file("trigrid_500.vtk");
    auto grid = UnstructuredGrid2D::vtk_read(gridfn);
    CavityFemSimpleWorker worker(grid, Re);
    worker.initialize_saver("cavity-fem-simple");
    worker.save_current_fields(0);

    // iterations loop
    size_t it = 0;
    for (it = 1; it < max_it; ++it) {
        double nrm = worker.init_step();
        std::cout << it - 1 << " " << nrm << std::endl;
        if (it > 0 && nrm < eps) {
            break;
        }
        worker.step();
        worker.save_current_fields(it);
    }
    CHECK(it == 24);
}
