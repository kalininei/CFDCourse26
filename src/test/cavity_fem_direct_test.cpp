#include "cfd/debug/printer.hpp"
#include "cfd/debug/saver.hpp"
#include "cfd/fem/elem2d/triangle_linear.hpp"
#include "cfd/fem/fem_assembler.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/umfpack_solver.hpp"
#include "cfd/numeric_integration/triangle_quadrature.hpp"
#include "cfd26_test.hpp"
#include "utils/filesystem.hpp"
#include "utils/vecmat.hpp"
#include <list>

using namespace cfd;

namespace {

struct CavityFemDirectWorker {
    CavityFemDirectWorker(const IGrid& grid, double Re);
    void initialize_saver(std::string stem);
    double init_step();
    void step();
    void save_current_fields(size_t iter);

private:
    const IGrid& grid_;
    const double Re_;
    FemAssembler fem_;

    CsrMatrix mat_;
    std::vector<double> rhs_;
    std::vector<double> uvp_;

    double* u_;
    double* v_;
    double* p_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_;

    static FemAssembler build_fem(const IGrid& grid);
};

CavityFemDirectWorker::CavityFemDirectWorker(const IGrid& grid, double Re)
    : grid_(grid),
      Re_(Re),
      fem_(build_fem(grid_)) {

    uvp_ = std::vector<double>(fem_.n_bases(), 0);
    rhs_ = std::vector<double>(fem_.n_bases(), 0);
    u_ = uvp_.data();
    v_ = u_ + grid_.n_points() + grid_.n_cells();
    p_ = v_ + grid_.n_points() + grid_.n_cells();
}

FemAssembler CavityFemDirectWorker::build_fem(const IGrid& grid) {
    constexpr int GAUSS_POWER = 3;
    const size_t iu0 = 0;
    const size_t iv0 = grid.n_points() + grid.n_cells();
    const size_t ip0 = 2 * (grid.n_points() + grid.n_cells());
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    // elements
    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        std::vector<size_t> ipoints = grid.tab_cell_point(icell);
        FemElement el;
        el.geometry = std::make_shared<TriangleLinearGeometry>(grid.point(ipoints[0]), grid.point(ipoints[1]),
                                                               grid.point(ipoints[2]));
        std::vector<std::shared_ptr<IElementBasis>> basis_list{std::make_shared<TriangleLinearBubbleBasis>(), // u
                                                               std::make_shared<TriangleLinearBubbleBasis>(), // v
                                                               std::make_shared<TriangleLinearBasis>()};      // p
        el.basis = std::make_shared<CompaundBasis>(basis_list);
        el.quadrature = quadrature_triangle_gauss<GAUSS_POWER>();
        elements.push_back(el);

        std::vector<size_t> elem_bases = {
            iu0 + ipoints[0], iu0 + ipoints[1], iu0 + ipoints[2], iu0 + grid.n_points() + icell, // u
            iv0 + ipoints[0], iv0 + ipoints[1], iv0 + ipoints[2], iv0 + grid.n_points() + icell, // v
            ip0 + ipoints[0], ip0 + ipoints[1], ip0 + ipoints[2]                                 // p
        };
        tab_elem_basis.push_back(elem_bases);
    }
    size_t n_bases = 3 * grid.n_points() + 2 * grid.n_cells();
    FemAssembler ret(n_bases, elements, tab_elem_basis);
    return ret;
}

void CavityFemDirectWorker::initialize_saver(std::string stem) {
    writer_.reset(new VtkUtils::TimeSeriesWriter(stem));
};

double CavityFemDirectWorker::init_step() {
    mat_.set_stencil(fem_.stencil());
    rhs_.resize(fem_.n_bases(), 0.0);

    for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
        const FemElement& el = fem_.element(ielem);
        std::vector<double> local_uvp_old = fem_.local_vector(ielem, uvp_);

        auto fun = [&](Point p) -> std::vector<double> {
            const size_t n = el.basis->size();
            std::vector<double> ret(n * n, 0.0);
            auto val = FemElementValue(&el, p);
            // convection velocity
            const Vector a = {val.subrange_interpolate(0, 4, local_uvp_old),
                              val.subrange_interpolate(4, 8, local_uvp_old)};

            // u
            for (size_t ibas = 0; ibas < 4; ++ibas) {
                Vector g1 = val.grad_phi(ibas);
                for (size_t jbas = 0; jbas < 4; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // -1/Re * laplace(u)
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] += (1.0 / Re_) * dot_product(g1, g2) * val.modj();
                    // a * grad(u)
                    ret[k1] += dot_product(a, g2) * val.phi(ibas) * val.modj();
                }
                for (size_t jbas = 8; jbas < 11; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // dp/dx
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] += g2.x * val.phi(ibas) * val.modj();
                }
            }
            // v
            for (size_t ibas = 4; ibas < 8; ++ibas) {
                Vector g1 = val.grad_phi(ibas);
                for (size_t jbas = 4; jbas < 8; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // -1/Re * laplace(v)
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] += (1.0 / Re_) * dot_product(g1, g2) * val.modj();
                    // a * grad(v)
                    ret[k1] += dot_product(a, g2) * val.phi(ibas) * val.modj();
                }
                for (size_t jbas = 8; jbas < 11; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    // dp/dy
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] += g2.y * val.phi(ibas) * val.modj();
                }
            }
            // p
            for (size_t ibas = 8; ibas < 11; ++ibas) {
                // du/dx
                for (size_t jbas = 0; jbas < 4; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] += g2.x * val.phi(ibas) * val.modj();
                }
                // dv/dx
                for (size_t jbas = 4; jbas < 8; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] += g2.y * val.phi(ibas) * val.modj();
                }
            }
            return ret;
        };
        auto local_mat = el.quadrature->integrate(fun);
        fem_.add_to_global_matrix(ielem, local_mat, mat_.vals());
    }

    // boundary conditions
    for (size_t i: grid_.boundary_points()) {
        // u = 0/1
        size_t u_idx = i;
        mat_.set_unit_row(u_idx);
        if (grid_.point(i).y > 0.99999) {
            rhs_[u_idx] = 1.0;
            uvp_[u_idx] = 1.0;
        } else {
            rhs_[u_idx] = 0.0;
            uvp_[u_idx] = 0.0;
        }
        // v = 0
        size_t v_idx = grid_.n_points() + grid_.n_cells() + i;
        mat_.set_unit_row(v_idx);
        rhs_[v_idx] = 0.0;
        uvp_[v_idx] = 0.0;
    }
    // p = 0 at ipoint=0
    mat_.set_unit_row(2 * grid_.n_points() + 2 * grid_.n_cells());
    rhs_[2 * grid_.n_points() + 2 * grid_.n_cells()] = 0.0;
    uvp_[2 * grid_.n_points() + 2 * grid_.n_cells()] = 0.0;

    return compute_residual(mat_, rhs_, uvp_);
}

void CavityFemDirectWorker::step() {
    UmfpackSolver slv;
    slv.set_matrix(mat_);
    slv.solve(rhs_, uvp_);
}

void CavityFemDirectWorker::save_current_fields(size_t iter) {
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

TEST_CASE("Cavity FEM-DIRECT", "[cavity-fem-direct]") {
    std::cout << std::endl << "--- cfd26_test [cavity-fem-direct] --- " << std::endl;

    // problem parameters
    double Re = 100;
    size_t max_it = 10'000;
    double eps = 1e-9;

    // worker initialization
    const std::string gridfn = test_directory_file("trigrid_500.vtk");
    auto grid = UnstructuredGrid2D::vtk_read(gridfn);
    CavityFemDirectWorker worker(grid, Re);
    worker.initialize_saver("cavity-fem-direct");
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
    CHECK(it == 12);
}
