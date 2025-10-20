#include "cfd/debug/printer.hpp"
#include "cfd/debug/saver.hpp"
#include "cfd/fem/elem1d/segment_linear.hpp"
#include "cfd/fem/elem2d/triangle_linear.hpp"
#include "cfd/fem/fem_assembler.hpp"
#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/regular_grid2d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd/numeric_integration/numeric_integration.hpp"
#include "cfd26_test.hpp"
#include "test/utils/filesystem.hpp"

using namespace cfd;

namespace {

///////////////////////////////////////////////////////////////////////////////
// ITestPoissonFemWorker
///////////////////////////////////////////////////////////////////////////////
struct ITestPoissonFemWorker {
    virtual ~ITestPoissonFemWorker() {}
    virtual double exact_solution(Point p) const = 0;
    virtual double exact_rhs(Point p) const = 0;
    virtual std::vector<size_t> dirichlet_bases() const {
        return grid_.boundary_points();
    }

    ITestPoissonFemWorker(const IGrid& grid, const FemAssembler& fem);
    virtual double solve();
    void save_vtk(const std::string& filename) const;

protected:
    const IGrid& grid_;
    FemAssembler fem_;
    std::vector<double> u_;

    virtual std::vector<double> element_load_vector(size_t ielem) const = 0;
    virtual std::vector<double> element_mass_matrix(size_t ielem) const = 0;
    virtual std::vector<double> element_stiffness_matrix(size_t ielem) const = 0;

    CsrMatrix approximate_lhs() const;
    std::vector<double> approximate_rhs() const;
    double compute_norm2() const;
};

ITestPoissonFemWorker::ITestPoissonFemWorker(const IGrid& grid, const FemAssembler& fem) : grid_(grid), fem_(fem) {}

double ITestPoissonFemWorker::solve() {
    // 1. build SLAE
    CsrMatrix mat = approximate_lhs();
    std::vector<double> rhs = approximate_rhs();
    // 2. solve SLAE
    AmgcMatrixSolver solver({{"precond.relax.type", "gauss_seidel"}});
    solver.set_matrix(mat);
    solver.solve(rhs, u_);
    // 3. compute norm2
    return compute_norm2();
}

void ITestPoissonFemWorker::save_vtk(const std::string& filename) const {
    // save grid
    grid_.save_vtk(filename);
    // save numerical solution
    VtkUtils::add_point_data(u_, "numerical", filename, grid_.n_points());
    // save exact solution
    std::vector<double> exact(grid_.n_points());
    for (size_t i = 0; i < grid_.n_points(); ++i) {
        exact[i] = exact_solution(grid_.point(i));
    }
    VtkUtils::add_point_data(exact, "exact", filename, grid_.n_points());
}

CsrMatrix ITestPoissonFemWorker::approximate_lhs() const {
    CsrMatrix ret(fem_.stencil());
    for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
        std::vector<double> local_stiff = element_stiffness_matrix(ielem);
        fem_.add_to_global_matrix(ielem, local_stiff, ret.vals());
    }
    // Dirichlet bc
    for (size_t ibas : dirichlet_bases()) {
        ret.set_unit_row(ibas);
    }
    return ret;
}

std::vector<double> ITestPoissonFemWorker::approximate_rhs() const {
    // mass matrix
    CsrMatrix mass(fem_.stencil());
    for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
        std::vector<double> local_mass = element_mass_matrix(ielem);
        fem_.add_to_global_matrix(ielem, local_mass, mass.vals());
    }
    // rhs = Mass * f
    std::vector<double> fvec(fem_.n_bases());
    for (size_t ibas = 0; ibas < fem_.n_bases(); ++ibas) {
        Point p = fem_.reference_point(ibas);
        fvec[ibas] = exact_rhs(p);
    }
    std::vector<double> ret = mass.mult_vec(fvec);
    // Dirichlet bc
    for (size_t ibas : dirichlet_bases()) {
        Point p = fem_.reference_point(ibas);
        ret[ibas] = exact_solution(p);
    }

    return ret;
}

double ITestPoissonFemWorker::compute_norm2() const {
    std::vector<double> load_vec(fem_.n_bases(), 0);
    for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
        std::vector<double> local_load = element_load_vector(ielem);
        fem_.add_to_global_vector(ielem, local_load, load_vec);
    }

    double integral = 0;
    double full_area = 0;
    for (size_t ibas = 0; ibas < fem_.n_bases(); ++ibas) {
        Point p = fem_.reference_point(ibas);
        double diff = u_[ibas] - exact_solution(p);
        integral += load_vec[ibas] * (diff * diff);
        full_area += load_vec[ibas];
    }
    return std::sqrt(integral / full_area);
}

////////////////////////////////////////////////////////////////////////////////
// ITestPoisson1FemWorker
////////////////////////////////////////////////////////////////////////////////
struct ITestPoisson1FemWorker : public ITestPoissonFemWorker {
    double exact_solution(Point p) const override {
        double x = p.x;
        return sin(10 * x * x);
    }
    double exact_rhs(Point p) const override {
        double x = p.x;
        return 400 * x * x * sin(10 * x * x) - 20 * cos(10 * x * x);
    }
    ITestPoisson1FemWorker(const IGrid& grid, const FemAssembler& fem) : ITestPoissonFemWorker(grid, fem) {}
};

///////////////////////////////////////////////////////////////////////////////
// TestPoissonLinearSegmentWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonLinearSegmentWorker : public ITestPoisson1FemWorker {
    TestPoissonLinearSegmentWorker(const IGrid& grid) : ITestPoisson1FemWorker(grid, build_fem(grid)) {}
    static FemAssembler build_fem(const IGrid& grid);

    std::vector<double> element_load_vector(size_t ielem) const override;
    std::vector<double> element_mass_matrix(size_t ielem) const override;
    std::vector<double> element_stiffness_matrix(size_t ielem) const override;
};

FemAssembler TestPoissonLinearSegmentWorker::build_fem(const IGrid& grid) {
    size_t n_bases = grid.n_points();
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    // elements
    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        std::vector<size_t> ipoints = grid.tab_cell_point(icell);
        Point p0 = grid.point(ipoints[0]);
        Point p1 = grid.point(ipoints[1]);

        auto geom = std::make_shared<SegmentLinearGeometry>(p0, p1);
        auto basis = std::make_shared<SegmentLinearBasis>();
        std::shared_ptr<const Quadrature> quad = nullptr;
        FemElement elem{geom, basis, quad};

        elements.push_back(elem);
        tab_elem_basis.push_back(ipoints);
    }

    return FemAssembler(n_bases, elements, tab_elem_basis);
}

std::vector<double> TestPoissonLinearSegmentWorker::element_load_vector(size_t ielem) const {
    double s = grid_.cell_volume(ielem) / 2.0;
    return {s, s};
}
std::vector<double> TestPoissonLinearSegmentWorker::element_mass_matrix(size_t ielem) const {
    double s = grid_.cell_volume(ielem) / 6;
    return {2 * s, s, s, 2 * s};
}
std::vector<double> TestPoissonLinearSegmentWorker::element_stiffness_matrix(size_t ielem) const {
    double s = 1.0 / grid_.cell_volume(ielem);
    return {s, -s, -s, s};
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// [poisson1-fem-linsegm]
///////////////////////////////////////////////////////////////////////////////
TEST_CASE("Poisson-fem 1D solver, linear segment elements", "[poisson1-fem-linsegm]") {
    std::cout << std::endl << "--- [poisson1-fem-linsegm] --- " << std::endl;
    Grid1D grid(0, 1, 10);
    TestPoissonLinearSegmentWorker worker(grid);
    double nrm = worker.solve();
    worker.save_vtk("poisson1_fem.vtk");
    std::cout << grid.n_cells() << " " << nrm << std::endl;
    CHECK(nrm == Approx(0.138156).margin(1e-6));
}
