#include "cfd/fvm/fvm_dfdn.hpp"
#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/fvm/fvm_gradient.hpp"
#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"
#include "utils/filesystem.hpp"

using namespace cfd;

namespace {

class ITestPoissonFvmWorker {
public:
    virtual ~ITestPoissonFvmWorker() = default;
    virtual double exact_solution(Point p) const = 0;
    virtual double exact_rhs(Point p) const = 0;

    // returns norm2(u - u_exact)
    double solve() {
        // 1. build SLAE
        CsrMatrix mat = approximate_lhs();
        std::vector<double> rhs = approximate_rhs();

        // 2. solve SLAE
        AmgcMatrixSolver solver;
        solver.set_matrix(mat);
        solver.solve(rhs, u_);

        // 3. compute norm2
        return compute_norm2();
    }

    // saves numerical and exact solution into the vtk format
    void save_vtk(const std::string& filename) {
        // save grid
        grid_->save_vtk(filename);

        // save numerical solution
        VtkUtils::add_cell_data(u_, "numerical", filename);

        // save exact solution
        std::vector<double> exact(grid_->n_cells());
        for (size_t i = 0; i < grid_->n_cells(); ++i) {
            exact[i] = exact_solution(grid_->cell_center(i));
        }
        VtkUtils::add_cell_data(exact, "exact", filename);
    }

    std::shared_ptr<const IGrid> grid() {
        return grid_;
    }

protected:
    ITestPoissonFvmWorker(std::shared_ptr<IGrid> grid) : grid_(grid){};

    std::shared_ptr<IGrid> grid_;
    std::vector<size_t> internal_faces_;
    struct DirichletFace {
        size_t iface;
        size_t icell;
        double value;
        Vector outer_normal;
    };
    std::vector<DirichletFace> dirichlet_faces_;
    std::vector<double> u_;

    void initialize() {
        // assemble face lists
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            size_t icell_negative = grid_->tab_face_cell(iface)[0];
            size_t icell_positive = grid_->tab_face_cell(iface)[1];
            if (icell_positive != INVALID_INDEX && icell_negative != INVALID_INDEX) {
                // internal faces list
                internal_faces_.push_back(iface);
            } else {
                // dirichlet faces list
                DirichletFace dir_face;
                dir_face.iface = iface;
                dir_face.value = exact_solution(grid_->face_center(iface));
                if (icell_positive == INVALID_INDEX) {
                    dir_face.icell = icell_negative;
                    dir_face.outer_normal = grid_->face_normal(iface);
                } else {
                    dir_face.icell = icell_positive;
                    dir_face.outer_normal = -grid_->face_normal(iface);
                }
                dirichlet_faces_.push_back(dir_face);
            }
        }
    }

    virtual CsrMatrix approximate_lhs() const {
        LodMatrix mat(grid_->n_cells());

        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            for (size_t iface: grid_->tab_cell_face(icell)) {
                std::array<size_t, 2> ij = grid_->tab_face_cell(iface);
                size_t jcell = (ij[0] == icell) ? ij[1] : ij[0];

                if (jcell != INVALID_INDEX) {
                    Point ci = grid_->cell_center(icell);
                    Point cj = grid_->cell_center(jcell);
                    double h = vector_abs(cj - ci);
                    double coef = grid_->face_area(iface) / h;
                    mat.add_value(icell, icell, coef);
                    mat.add_value(icell, jcell, -coef);
                } else {
                    Point gs = grid_->face_center(iface);
                    Point ci = grid_->cell_center(icell);
                    double h = vector_abs(gs - ci);
                    double coef = grid_->face_area(iface) / h;
                    mat.add_value(icell, icell, coef);
                }
            }
        }

        return mat.to_csr();
    }

    virtual std::vector<double> approximate_rhs() const {
        std::vector<double> rhs(grid_->n_cells(), 0.0);
        // internal
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double value = exact_rhs(grid_->cell_center(icell));
            double volume = grid_->cell_volume(icell);
            rhs[icell] = value * volume;
        }
        // dirichlet faces
        for (const DirichletFace& dir_face: dirichlet_faces_) {
            size_t icell = dir_face.icell;
            size_t iface = dir_face.iface;
            Point gs = grid_->face_center(iface);
            Point ci = grid_->cell_center(icell);
            double h = vector_abs(gs - ci);
            double coef = grid_->face_area(iface) / h;
            rhs[icell] += dir_face.value * coef;
        }
        return rhs;
    }

    double compute_norm2() const {
        double norm2 = 0;
        double full_area = 0;
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double diff = u_[icell] - exact_solution(grid_->cell_center(icell));
            norm2 += grid_->cell_volume(icell) * diff * diff;
            full_area += grid_->cell_volume(icell);
        }
        return std::sqrt(norm2 / full_area);
    }
};

/////////////////////////////////////////////////
// 1D Worker
/////////////////////////////////////////////////
class TestPoissonFvmWorker_1D : public ITestPoissonFvmWorker {
public:
    TestPoissonFvmWorker_1D(size_t n_cells) : ITestPoissonFvmWorker(std::make_shared<Grid1D>(0, 1, n_cells)) {
        initialize();
    }

    // u(x) = sin(10*x^2)
    double exact_solution(Point p) const override {
        double x = p.x;
        return sin(10 * x * x);
    }

    // -d^2 u(x)/d x^2
    double exact_rhs(Point p) const override {
        double x = p.x;
        return 400 * x * x * sin(10 * x * x) - 20 * cos(10 * x * x);
    }
};

} // namespace

TEST_CASE("Poisson 1D solver, Finite Volume Method", "[poisson1-fvm]") {
    std::cout << std::endl << "--- [poisson1-fvm] --- " << std::endl;

    // precalculated norm2 results for some n_cells values
    // used for CHECK procedures
    std::map<size_t, double> norm2_for_compare{
        {10, 0.106539},
        {100, 0.00101714},
        {1000, 1.01641e-05},
    };

    // loop over n_cells value
    for (size_t n_cells: {10, 20, 50, 100, 200, 500, 1000}) {
        // build test solver
        TestPoissonFvmWorker_1D worker(n_cells);

        // solve and find norm2
        double n2 = worker.solve();

        // save into poisson1_ncells={n_cells}.vtk
        worker.save_vtk("poisson1_fvm_n=" + std::to_string(n_cells) + ".vtk");

        // print (N_CELLS, NORM2) table entry
        std::cout << n_cells << " " << n2 << std::endl;

        // CHECK if result for this n_cells
        // presents in the norm2_for_compare dictionary
        auto found = norm2_for_compare.find(n_cells);
        if (found != norm2_for_compare.end()) {
            CHECK(n2 == Approx(found->second).margin(1e-6));
        }
    }
}

namespace {

/////////////////////////////////////////////////
// 2D Worker
/////////////////////////////////////////////////

class TestPoissonFvmWorker_2D : public ITestPoissonFvmWorker {
public:
    TestPoissonFvmWorker_2D(std::string grid_fn)
        : ITestPoissonFvmWorker(std::make_shared<UnstructuredGrid2D>(UnstructuredGrid2D::vtk_read(grid_fn, true))),
          ecol_(*grid()), grad_computer_(std::make_shared<LeastSquaresFvmFaceGradient>(*grid(), ecol_)) {
        initialize();
        u_.resize(ecol_.size());
        std::fill(u_.begin(), u_.end(), 0.0);
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

    CsrMatrix approximate_lhs() const override {
        LodMatrix mat(ecol_.size());
        // internal
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [i, j] = ecol_.tab_face_colloc(iface);
            double area = grid_->face_area(iface);
            double hij = vector_abs(ecol_.points[i] - ecol_.points[j]);
            double coef = area / hij;

            mat.add_value(i, i, coef);
            mat.add_value(i, j, -coef);
            mat.add_value(j, j, coef);
            mat.add_value(j, i, -coef);
        }
        // dirichlet faces
        for (auto& dir: dirichlet_faces_) {
            size_t icolloc = ecol_.boundary_colloc(dir.iface);
            mat.set_unit_row(icolloc);
        }

        return mat.to_csr();
    }

    std::vector<double> approximate_rhs() const override {
        std::vector<double> rhs(ecol_.size(), 0.0);

        // non-orhtogonality
        std::vector<Vector> gradu = grad_computer_->compute(u_);
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [i, j] = ecol_.tab_face_colloc(iface);
            double area = grid_->face_area(iface);
            Vector cij = ecol_.points[j] - ecol_.points[i];
            double hij = vector_abs(cij);
            Vector a = grid_->face_normal(iface) - cij / hij;
            double val = dot_product(gradu[iface], a) * area;

            rhs[i] += val;
            rhs[j] -= val;
        }
        // internal
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double value = exact_rhs(grid_->cell_center(icell));
            double volume = grid_->cell_volume(icell);
            rhs[icell] += value * volume;
        }
        // dirichlet faces
        for (const DirichletFace& dir: dirichlet_faces_) {
            size_t icolloc = ecol_.boundary_colloc(dir.iface);
            rhs[icolloc] = dir.value;
        }

        return rhs;
    }

private:
    FvmExtendedCollocations ecol_;
    std::shared_ptr<IFvmFaceGradient> grad_computer_;
};

} // namespace

TEST_CASE("Poisson 2D solver, unstructured", "[poisson2-fvm]") {
    std::cout << std::endl << "--- [poisson2-fvm] --- " << std::endl;

    TestPoissonFvmWorker_2D worker(test_directory_file("tetragrid_500.vtk"));
    double n2;
    for (size_t i = 0; i < 10; ++i) {
        n2 = worker.solve();
        std::cout << n2 << std::endl;
    }
    CHECK(n2 == Approx(0.0411111).margin(1e-6));
}

namespace {

class RadialGrid1D : public Grid1D {
public:
    RadialGrid1D(double r0, double r1, size_t n_cells) : Grid1D(r0, r1, n_cells) {}

    double face_area(size_t iface) const override {
        return 2 * M_PI * face_center(iface).x * Grid1D::face_area(iface);
    }
    double cell_volume(size_t icell) const override {
        return 2 * M_PI * cell_center(icell).x * Grid1D::cell_volume(icell);
    }
};

/////////////////////////////////////////////////
// Radial1D Worker
/////////////////////////////////////////////////

class TestPoissonFvmWorker_Radial1D : public ITestPoissonFvmWorker {
public:
    TestPoissonFvmWorker_Radial1D(double r0, double lambda0, size_t n_cells)
        : ITestPoissonFvmWorker(std::make_shared<RadialGrid1D>(r0, 1.0 + r0, n_cells)), ecol_(*grid()), r0_(r0) {

        // lambda
        lambda_.resize(ecol_.size());
        for (size_t i = 0; i < ecol_.size(); ++i) {
            lambda_[i] = (ecol_.points[i].x < r0 + 0.5) ? lambda0 : 1.0;
        }

        // solution vector
        u_.resize(ecol_.size());
        std::fill(u_.begin(), u_.end(), 0.0);

        // initialize
        initialize();
    }
    double exact_solution(Point p) const override {
        const double r = p.x;
        const double r0 = r0_;
        const double r1 = r0_ + 1.0;
        const double rk = r0_ + 0.5;
        const double k0 = lambda_.front();
        const double k1 = lambda_.back();
        const double a = 1.0;
        const double b = 0.0;
        const double B = std::log(rk / r0) / k0 + std::log(r1 / rk) / k1;
        if (r < rk) {
            return a + (b - a) / (B * k0) * std::log(r / r0);
        } else {
            return b + (b - a) / (B * k1) * std::log(r / r1);
        }
    }
    double exact_rhs(Point) const override {
        return 0;
    }

    CsrMatrix approximate_lhs() const override {
        LodMatrix mat(ecol_.size());
        // internal
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [i, j] = ecol_.tab_face_colloc(iface);
            double area = grid_->face_area(iface);
            double hij = vector_abs(ecol_.points[i] - ecol_.points[j]);
            double coef = face_lambda(iface) * area / hij;

            mat.add_value(i, i, coef);
            mat.add_value(i, j, -coef);
            mat.add_value(j, j, coef);
            mat.add_value(j, i, -coef);
        }
        // dirichlet faces
        for (auto& dir: dirichlet_faces_) {
            size_t icolloc = ecol_.boundary_colloc(dir.iface);
            mat.set_unit_row(icolloc);
        }

        return mat.to_csr();
    }

    std::vector<double> approximate_rhs() const override {
        std::vector<double> rhs(ecol_.size(), 0.0);
        // internal
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double value = exact_rhs(grid_->cell_center(icell));
            double volume = grid_->cell_volume(icell);
            rhs[icell] += value * volume;
        }
        // dirichlet faces
        for (const DirichletFace& dir: dirichlet_faces_) {
            size_t icolloc = ecol_.boundary_colloc(dir.iface);
            rhs[icolloc] = dir.value;
        }

        return rhs;
    }

    double face_lambda(size_t iface) const {
        auto [i, j] = ecol_.tab_face_colloc(iface);
        return (lambda_[i] + lambda_[j]) / 2.0;
    }

private:
    FvmExtendedCollocations ecol_;
    const double r0_;
    std::vector<double> lambda_;
};

} // namespace

TEST_CASE("Poisson radial solver 1D", "[poisson1-fvm-radial]") {
    std::cout << std::endl << "--- [poisson1-fvm-radial] --- " << std::endl;

    size_t n_cells = 10;
    TestPoissonFvmWorker_Radial1D worker(0.05, 10, n_cells);
    double n2 = worker.solve();
    std::cout << n_cells << " " << n2 << std::endl;
    CHECK(n2 == Approx(0.0252116).margin(1e-6));
}
