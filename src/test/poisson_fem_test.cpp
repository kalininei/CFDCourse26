#include "cfd/debug/printer.hpp"
#include "cfd/debug/saver.hpp"
#include "cfd/fem/elem1d/segment_cubic.hpp"
#include "cfd/fem/elem1d/segment_linear.hpp"
#include "cfd/fem/elem1d/segment_quadratic.hpp"
#include "cfd/fem/elem2d/quadrangle_linear.hpp"
#include "cfd/fem/elem2d/quadrangle_quadratic.hpp"
#include "cfd/fem/elem2d/triangle_cubic.hpp"
#include "cfd/fem/elem2d/triangle_linear.hpp"
#include "cfd/fem/elem2d/triangle_quadratic.hpp"
#include "cfd/fem/fem_assembler.hpp"
#include "cfd/fem/fem_boundary.hpp"
#include "cfd/fem/fem_sorted_cell_info.hpp"
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
    virtual std::vector<size_t> neumann_faces() const {
        return {};
    };

    ITestPoissonFemWorker(const IGrid& grid, const FemAssembler& fem);
    virtual double solve();
    void save_vtk(const std::string& filename) const;

    const std::vector<double>& u() const {
        return u_;
    }

protected:
    const IGrid& grid_;
    FemAssembler fem_;
    std::vector<double> u_;

    virtual std::vector<double> element_load_vector(size_t ielem) const;
    virtual std::vector<double> element_mass_matrix(size_t ielem) const;
    virtual std::vector<double> element_stiffness_matrix(size_t ielem) const;
    virtual std::vector<double> face_neumann_vector(size_t) const {
        _THROW_NOT_IMP_;
    }

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
    AmgcMatrixSolver solver({{"precond.relax.type", "gauss_seidel"}, {"solver.tol", "1e-12"}});
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
    for (size_t ibas: dirichlet_bases()) {
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
    // Neumann bc
    for (size_t iface: neumann_faces()) {
        std::vector<double> local_neu = face_neumann_vector(iface);
        fem_.boundary_add_to_global_vector(iface, local_neu, ret);
    }
    // Dirichlet bc
    for (size_t ibas: dirichlet_bases()) {
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

std::vector<double> ITestPoissonFemWorker::element_load_vector(size_t ielem) const {
    const FemElement& el = fem_.element(ielem);

    auto fun = [&el](Point p) -> std::vector<double> {
        std::vector<double> ret(el.basis->size(), 0.0);
        auto val = FemElementValue(&el, p);
        for (size_t ibas = 0; ibas < ret.size(); ++ibas) {
            ret[ibas] = val.phi(ibas) * val.modj();
        }
        return ret;
    };

    return el.quadrature->integrate(fun);
}

std::vector<double> ITestPoissonFemWorker::element_mass_matrix(size_t ielem) const {
    const FemElement& el = fem_.element(ielem);

    auto fun = [&el](Point p) -> std::vector<double> {
        const size_t n = el.basis->size();
        std::vector<double> ret(n * n, 0.0);
        auto val = FemElementValue(&el, p);

        for (size_t ibas = 0; ibas < n; ++ibas) {
            for (size_t jbas = 0; jbas <= ibas; ++jbas) {
                const size_t k1 = ibas * n + jbas;
                ret[k1] = val.phi(ibas) * val.phi(jbas) * val.modj();
                if (ibas != jbas) {
                    const size_t k2 = jbas * n + ibas;
                    ret[k2] = ret[k1];
                }
            }
        }

        return ret;
    };

    return el.quadrature->integrate(fun);
}

std::vector<double> ITestPoissonFemWorker::element_stiffness_matrix(size_t ielem) const {
    const FemElement& el = fem_.element(ielem);

    auto fun = [&el](Point p) -> std::vector<double> {
        const size_t n = el.basis->size();
        std::vector<double> ret(n * n, 0.0);
        auto val = FemElementValue(&el, p);

        for (size_t ibas = 0; ibas < n; ++ibas) {
            for (size_t jbas = 0; jbas <= ibas; ++jbas) {
                const size_t k1 = ibas * n + jbas;
                Vector g1 = val.grad_phi(ibas);
                Vector g2 = val.grad_phi(jbas);
                ret[k1] = dot_product(g1, g2) * val.modj();
                if (ibas != jbas) {
                    const size_t k2 = jbas * n + ibas;
                    ret[k2] = ret[k1];
                }
            }
        }

        return ret;
    };

    return el.quadrature->integrate(fun);
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

namespace {
///////////////////////////////////////////////////////////////////////////////
// TestPoissonQuadraticTriangleWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonQuadraticTriangleWorker : public ITestPoissonFemWorker {

    TestPoissonQuadraticTriangleWorker(const IGrid& grid) : ITestPoissonFemWorker(grid, build_fem(grid)) {}
    static FemAssembler build_fem(const IGrid& grid);

    std::vector<size_t> dirichlet_bases() const override {
        std::vector<size_t> ret = grid_.boundary_points();
        for (size_t iface: grid_.boundary_faces()) {
            ret.push_back(grid_.n_points() + iface);
        }
        return ret;
    }

    double exact_solution(Point p) const override {
        double x = p.x;
        double y = p.y;
        return cos(10 * x * x) * sin(10 * y) + sin(10 * x * x) * cos(10 * x);
    }

    double exact_rhs(Point p) const override {
        double x = p.x;
        double y = p.y;
        return (20 * sin(10 * x * x) + (400 * x * x + 100) * cos(10 * x * x)) * sin(10 * y) +
               (400 * x * x + 100) * cos(10 * x) * sin(10 * x * x) +
               (400 * x * sin(10 * x) - 20 * cos(10 * x)) * cos(10 * x * x);
    }
};

FemAssembler TestPoissonQuadraticTriangleWorker::build_fem(const IGrid& grid) {
    size_t n_bases = grid.n_points() + grid.n_faces();
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    // elements
    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        if (grid.tab_cell_point(icell).size() != 3) {
            // only triangle grid is allowed
            _THROW_NOT_IMP_;
        }
        // this object provides matched sorted points and faces list for the given cell
        // so that info.ifaces[0] lies between info.ipoints[0] and info.ipoints[1]
        PolygonElementInfo info(grid, icell);

        Point p0 = grid.point(info.ipoints[0]);
        Point p1 = grid.point(info.ipoints[1]);
        Point p2 = grid.point(info.ipoints[2]);

        auto geom = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
        auto basis = std::make_shared<TriangleQuadraticBasis>();
        auto quad = quadrature_triangle_gauss6();
        FemElement elem{geom, basis, quad};

        elements.push_back(elem);
        std::vector<size_t> tab = info.ipoints;
        for (size_t i=0; i<info.n_points(); ++i){
            size_t iface = info.ifaces[i];
            tab.push_back(grid.n_points() + iface);
        }
        tab_elem_basis.push_back(tab);
    }

    return FemAssembler(n_bases, elements, tab_elem_basis);
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// [poisson2-fem-quadtri]
///////////////////////////////////////////////////////////////////////////////
TEST_CASE("Poisson-fem 2D solver, quadratic triangle elements", "[poisson2-fem-quadtri]") {
    std::cout << std::endl << "--- [poisson2-fem-quadtri] --- " << std::endl;

    const std::string gridfn = test_directory_file("trigrid_500.vtk");
    auto grid = UnstructuredGrid2D::vtk_read(gridfn, true);
    TestPoissonQuadraticTriangleWorker worker(grid);
    double nrm = worker.solve();
    worker.save_vtk("poisson2.vtk");
    std::cout << grid.n_cells() << " " << nrm << std::endl;
    CHECK(nrm == Approx(0.0048626614).margin(1e-6));
}

namespace {

///////////////////////////////////////////////////////////////////////////////
// TestPoissonRadialWorker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonRadialWorker : public ITestPoissonFemWorker {

    TestPoissonRadialWorker(const Grid1D& grid)
        : ITestPoissonFemWorker(grid, build_fem(grid)), r0_(grid.point(0).x), r1_(grid.point(grid.n_points() - 1).x),
          q_(-1.0) {}

    static FemAssembler build_fem(const IGrid& grid);

    double exact_solution(Point p) const override {
        double r = p.x;
        return q_ * r0_ * std::log(r / r1_);
    }

    double exact_rhs(Point) const override {
        return 0;
    }

    std::vector<size_t> dirichlet_bases() const override {
        return {grid_.n_points() - 1};
    }

    std::vector<size_t> neumann_faces() const override {
        return {0};
    };

    std::vector<double> element_load_vector(size_t ielem) const override;
    std::vector<double> element_mass_matrix(size_t ielem) const override;
    std::vector<double> element_stiffness_matrix(size_t ielem) const override;
    std::vector<double> face_neumann_vector(size_t iface) const override;

private:
    const double r0_;
    const double r1_;
    const double q_;
};

FemAssembler TestPoissonRadialWorker::build_fem(const IGrid& grid) {
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
        auto quad = quadrature_segment_gauss6();
        FemElement elem{geom, basis, quad};

        elements.push_back(elem);

        std::vector<size_t> tab{ipoints[0], ipoints[1]};
        tab_elem_basis.push_back(tab);
    }
    FemAssembler ret(n_bases, elements, tab_elem_basis);

    // boundary elements
    ret.add_boundary_element(0, 0, build_boundary_element(elements[0], nullptr, [](Point) { return Point{-1}; }));

    return ret;
}

std::vector<double> TestPoissonRadialWorker::element_load_vector(size_t ielem) const {
    const FemElement& el = fem_.element(ielem);

    auto fun = [&el](Point p) -> std::vector<double> {
        std::vector<double> ret(el.basis->size(), 0.0);
        auto val = FemElementValue(&el, p);
        double r = el.geometry->to_physical(p).x;

        for (size_t ibas = 0; ibas < ret.size(); ++ibas) {
            ret[ibas] = 2 * M_PI * val.phi(ibas) * r * val.modj();
        }
        return ret;
    };

    return el.quadrature->integrate(fun);
}

std::vector<double> TestPoissonRadialWorker::element_mass_matrix(size_t ielem) const {
    const FemElement& el = fem_.element(ielem);

    auto fun = [&el](Point p) -> std::vector<double> {
        const size_t n = el.basis->size();
        std::vector<double> ret(n * n, 0.0);
        auto val = FemElementValue(&el, p);
        double r = el.geometry->to_physical(p).x;

        for (size_t ibas = 0; ibas < n; ++ibas) {
            for (size_t jbas = 0; jbas <= ibas; ++jbas) {
                const size_t k1 = ibas * n + jbas;
                ret[k1] = 2 * M_PI * val.phi(ibas) * val.phi(jbas) * r * val.modj();
                if (ibas != jbas) {
                    const size_t k2 = jbas * n + ibas;
                    ret[k2] = ret[k1];
                }
            }
        }

        return ret;
    };

    return el.quadrature->integrate(fun);
}

std::vector<double> TestPoissonRadialWorker::element_stiffness_matrix(size_t ielem) const {
    const FemElement& el = fem_.element(ielem);

    auto fun = [&el](Point p) -> std::vector<double> {
        const size_t n = el.basis->size();
        std::vector<double> ret(n * n, 0.0);
        auto val = FemElementValue(&el, p);
        double r = el.geometry->to_physical(p).x;

        for (size_t ibas = 0; ibas < n; ++ibas) {
            for (size_t jbas = 0; jbas <= ibas; ++jbas) {
                const size_t k1 = ibas * n + jbas;
                Vector g1 = val.grad_phi(ibas);
                Vector g2 = val.grad_phi(jbas);
                ret[k1] = 2 * M_PI * dot_product(g1, g2) * r * val.modj();
                if (ibas != jbas) {
                    const size_t k2 = jbas * n + ibas;
                    ret[k2] = ret[k1];
                }
            }
        }

        return ret;
    };

    return el.quadrature->integrate(fun);
}

std::vector<double> TestPoissonRadialWorker::face_neumann_vector(size_t iface) const {
    const FemElement& el = fem_.boundary_element(iface);
    const size_t n = el.basis->size();
    std::vector<double> ret(n);

    auto val = FemElementValue(&el, {});
    for (size_t ibas = 0; ibas < n; ++ibas) {
        ret[ibas] = -2 * M_PI * q_ * r0_ * val.phi(ibas) * val.modj();
    }

    return ret;
};

} // namespace

///////////////////////////////////////////////////////////////////////////////
// [poisson1-fem-linsegm-radial]
///////////////////////////////////////////////////////////////////////////////
TEST_CASE("Poisson-fem 1D radial solver, linear segment elements", "[poisson1-fem-radial]") {
    std::cout << std::endl << "--- [poisson1-fem-radial] --- " << std::endl;
    const double r0 = 0.05;
    const double r1 = 1.05;
    const size_t n = 10;
    Grid1D grid(r0, r1, n);
    TestPoissonRadialWorker worker(grid);
    double nrm = worker.solve();
    worker.save_vtk("poisson1_fem_radial.vtk");
    double d = std::abs(worker.exact_solution({r0}) - worker.u()[0]);
    std::cout << grid.n_cells() << " " << nrm << " " << d << std::endl;
    CHECK(nrm == Approx(0.000528306).margin(1e-6));
}

namespace {

///////////////////////////////////////////////////////////////////////////////
// TestPoissonRadial2Worker
///////////////////////////////////////////////////////////////////////////////
struct TestPoissonRadial2Worker : public ITestPoissonFemWorker {

    TestPoissonRadial2Worker(const IGrid& grid, double r0, double r1)
        : ITestPoissonFemWorker(grid, build_fem(grid, r0)), r0_(r0), r1_(r1), q_(-1.0) {}

    static FemAssembler build_fem(const IGrid& grid, double r0);

    double exact_solution(Point p) const override {
        double r = std::sqrt(p.x * p.x + p.y * p.y);
        return q_ * r0_ * std::log(r / r1_);
    }

    double exact_rhs(Point) const override {
        return 0;
    }

    std::vector<size_t> dirichlet_bases() const override {
        std::set<size_t> ret;
        for (size_t iface: grid_.boundary_faces()) {
            if (vector_abs(grid_.face_center(iface)) > (r0_ + r1_) / 2.0) {
                // iface is a face from the outer boundary
                for (size_t ipoint: grid_.tab_face_point(iface)) {
                    ret.insert(ipoint);
                }
            }
        }
        return std::vector<size_t>(ret.begin(), ret.end());
    }

    std::vector<size_t> neumann_faces() const override {
        std::vector<size_t> ret;
        for (size_t iface: grid_.boundary_faces()) {
            if (vector_abs(grid_.face_center(iface)) <= r0_) {
                ret.push_back(iface);
            }
        }
        return ret;
    };

    std::vector<double> face_neumann_vector(size_t iface) const override;

private:
    const double r0_;
    const double r1_;
    const double q_;
};

FemAssembler TestPoissonRadial2Worker::build_fem(const IGrid& grid, double r0) {
    constexpr int GAUSS_POWER = 6;
    constexpr double BND_EPS = 1e-8;
    size_t n_bases = grid.n_points();
    std::vector<FemElement> elements;
    std::vector<std::vector<size_t>> tab_elem_basis;

    // elements
    std::map<size_t, size_t> bnd_face_cell;
    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        FemElement el;
        PolygonElementInfo info(grid, icell);
        info.start_from_boundary(); // if this is a boundary cell, then info.iface[0] is a boundary face

        Point p0 = grid.point(info.ipoints[0]);
        Point p1 = grid.point(info.ipoints[1]);
        Point p2 = grid.point(info.ipoints[2]);
        if (info.n_points() == 3) {
            // P-element for 3 points cell. Guaranteed to be not inner boundary cell.
            el.geometry = std::make_shared<TriangleLinearGeometry>(p0, p1, p2);
            el.basis = std::make_shared<TriangleLinearBasis>();
            el.quadrature = quadrature_triangle_gauss<GAUSS_POWER>();
        } else {
            // Q-element for 4 points cell
            Point p3 = grid.point(info.ipoints[3]);
            // if inner boundary cell
            if (vector_abs(grid.face_center(info.ifaces[0])) < r0 + BND_EPS) {
                // geometry transformation shape based on {(-1, -1), (1, -1), (1, 1), (-1, 1), (0, -1)}
                auto s_basis = std::make_shared<Quadrangle_Sqr0_Linear123>();
                // bubble point on the arc between p0 and p1
                Point p4 = r0 * vector_normalize(p0 + p1);
                el.geometry = build_geometry_from_basis(s_basis, {p0, p1, p2, p3, p4});
            } else {
                // otherwise use linear transformation
                el.geometry = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
            }
            el.basis = std::make_shared<QuadrangleLinearBasis>();
            el.quadrature = quadrature_square_gauss<GAUSS_POWER>();
        }
        elements.push_back(el);
        std::vector<size_t> tab = info.ipoints;
        tab_elem_basis.push_back(tab);
    }
    FemAssembler ret(n_bases, elements, tab_elem_basis);

    // create boundary element on the inner boundary
    for (auto iface: grid.boundary_faces()) {
        // if inner
        if (vector_abs(grid.face_center(iface)) < r0 + BND_EPS) {
            // get cell index
            auto [icell, c2] = grid.tab_face_cell(iface);
            if (icell == INVALID_INDEX) {
                icell = c2;
            }
            // assemble boundary element
            auto belem = build_boundary_element(elements[icell], quadrature_segment_gauss<GAUSS_POWER>(),
                                                // xi(tau) = tau, eta(tau) = -1
                                                [](Point xi) { return Point{xi.x, -1.0, 0.0}; },
                                                // Jtau = (dxi/dtau, 0, 0)
                                                {1.0, 0.0, 0.0});
            ret.add_boundary_element(iface, icell, belem);
        }
    }

    return ret;
}

std::vector<double> TestPoissonRadial2Worker::face_neumann_vector(size_t iface) const {
    const FemElement& el = fem_.boundary_element(iface);

    auto fun = [&](Point p) -> std::vector<double> {
        const size_t n = el.basis->size();
        std::vector<double> ret(n, 0.0);
        auto val = FemElementValue(&el, p);
        for (size_t ibas = 0; ibas < n; ++ibas) {
            ret[ibas] = -q_ * val.phi(ibas) * val.modj();
        }
        return ret;
    };

    return el.quadrature->integrate(fun);
};

} // namespace

///////////////////////////////////////////////////////////////////////////////
// [poisson2-fem-radial]
///////////////////////////////////////////////////////////////////////////////
TEST_CASE("Poisson-fem 2D radial solver, linear elements", "[poisson2-fem-radial]") {
    std::cout << std::endl << "--- [poisson2-fem-radial] --- " << std::endl;

    const double r0 = 0.05;
    const double r1 = 1.05;
    const std::string gridfn = test_directory_file("radial_500.vtk");
    auto grid = UnstructuredGrid2D::vtk_read(gridfn, true);
    TestPoissonRadial2Worker worker(grid, r0, r1);
    double nrm = worker.solve();
    worker.save_vtk("poisson2.vtk");
    std::cout << grid.n_cells() << " " << nrm << std::endl;
    CHECK(nrm == Approx(0.0003809977).margin(1e-6));
}
