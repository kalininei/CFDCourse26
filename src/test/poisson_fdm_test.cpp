#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Poisson 1D solver
///////////////////////////////////////////////////////////////////////////////

namespace {
class TestPoisson1Worker {
public:
    // u(x) = sin(10*x^2)
    double exact_solution(double x) const {
        return sin(10 * x * x);
    }

    // -d^2 u(x)/d x^2
    double exact_rhs(double x) const {
        return 400 * x * x * sin(10 * x * x) - 20 * cos(10 * x * x);
    }

    TestPoisson1Worker(size_t n_cells) : grid_(0, 1, n_cells) {}

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
        grid_.save_vtk(filename);

        // save numerical solution
        VtkUtils::add_point_data(u_, "numerical", filename);

        // save exact solution
        std::vector<double> exact(grid_.n_points());
        for (size_t i = 0; i < grid_.n_points(); ++i) {
            exact[i] = exact_solution(grid_.point(i).x);
        }
        VtkUtils::add_point_data(exact, "exact", filename);
    }

private:
    const Grid1D grid_;
    std::vector<double> u_;

    CsrMatrix approximate_lhs() const {
        // constant h = x[1] - x[0]
        double h = grid_.point(1).x - grid_.point(0).x;

        // fill using 'easy-to-construct' sparse matrix format
        LodMatrix mat(grid_.n_points());
        mat.add_value(0, 0, 1);
        mat.add_value(grid_.n_points() - 1, grid_.n_points() - 1, 1);
        double diag = 2.0 / h / h;
        double nondiag = -1.0 / h / h;
        for (size_t i = 1; i < grid_.n_points() - 1; ++i) {
            mat.add_value(i, i - 1, nondiag);
            mat.add_value(i, i + 1, nondiag);
            mat.add_value(i, i, diag);
        }

        // return 'easy-to-use' sparse matrix format
        return mat.to_csr();
    }

    std::vector<double> approximate_rhs() const {
        std::vector<double> ret(grid_.n_points());
        ret[0] = exact_solution(grid_.point(0).x);
        ret[grid_.n_points() - 1] = exact_solution(grid_.point(grid_.n_points() - 1).x);
        for (size_t i = 1; i < grid_.n_points() - 1; ++i) {
            ret[i] = exact_rhs(grid_.point(i).x);
        }
        return ret;
    }

    double compute_norm2() const {
        // weights
        double h = grid_.point(1).x - grid_.point(0).x;
        std::vector<double> w(grid_.n_points(), h);
        w[0] = w[grid_.n_points() - 1] = h / 2;

        // sum
        double sum = 0;
        for (size_t i = 0; i < grid_.n_points(); ++i) {
            double diff = u_[i] - exact_solution(grid_.point(i).x);
            sum += w[i] * diff * diff;
        }

        double len = grid_.point(grid_.n_points() - 1).x - grid_.point(0).x;
        return std::sqrt(sum / len);
    }
};
} // namespace

TEST_CASE("Poisson 1D solver, Finite Difference Method", "[poisson1-fdm]") {
    std::cout << std::endl << "--- [poisson1-fdm] --- " << std::endl;

    // precalculated norm2 results for some n_cells values
    // used for CHECK procedures
    std::map<size_t, double> norm2_for_compare{
        {10, 0.179124},
        {100, 0.00158055},
        {1000, 1.57849e-05},
    };

    // loop over n_cells value
    for (size_t n_cells : {10, 20, 50, 100, 200, 500, 1000}) {
        // build test solver
        TestPoisson1Worker worker(n_cells);

        // solve and find norm2
        double n2 = worker.solve();

        // save into poisson1_ncells={n_cells}.vtk
        worker.save_vtk("poisson1_fdm_n=" + std::to_string(n_cells) + ".vtk");

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
