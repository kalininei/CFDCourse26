#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"
#include <ranges>

using namespace cfd;

namespace {

class ATestTransport1Worker {
public:
    virtual ~ATestTransport1Worker() {}

    static double init_solution(double x) {
        // return (x >= -0.1 && x <= 0.1) ? 1.0 : 0.0;
        constexpr double sigma = 0.1;
        return exp(-x * x / sigma / sigma);
    }
    double exact_solution(double x) const {
        return init_solution(x - time_);
    }

    ATestTransport1Worker(size_t n_cells, double tau)
        : grid_(0, 1, n_cells), h_(1.0 / static_cast<double>(n_cells)), tau_(tau), u_(grid_.n_points()) {
        for (size_t i = 0; i < grid_.n_points(); ++i) {
            u_[i] = init_solution(grid_.point(i).x);
        }
    }

    void step() {
        time_ += tau_;
        impl_step();
    }

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

    double current_time() const {
        return time_;
    }

    double compute_norm2() {
        std::vector<double> diff(u_.size());
        for (size_t i = 0; i < grid_.n_points(); ++i) {
            diff[i] = u_[i] - exact_solution(grid_.point(i).x);
        }
        return compute_norm2(diff);
    }

    double compute_norm2(const std::vector<double>& v) {
        // weights
        std::vector<double> w(grid_.n_points(), h_);
        w.front() = w.back() = h_ / 2;

        // sum
        double sum = 0;
        for (size_t i = 0; i < w.size(); ++i) {
            sum += w[i] * v[i] * v[i];
        }

        double len = grid_.point(grid_.n_points() - 1).x - grid_.point(0).x;
        return std::sqrt(sum / len);
    }

protected:
    Grid1D grid_;
    const double h_;
    const double tau_;
    std::vector<double> u_;
    double time_ = 0;

private:
    virtual void impl_step() = 0;
};

enum class Limiter {
    UPWIND,
    CENTRAL,
    MINMOD,
};

template<Limiter L>
double limiter([[maybe_unused]] double r) {
    if constexpr (L == Limiter::UPWIND) {
        return 0.0;
    } else if constexpr (L == Limiter::CENTRAL) {
        return 1.0;
    } else if constexpr (L == Limiter::MINMOD) {
        return std::max(0.0, std::min(1.0, r));
    } else {
        static_assert(false);
    }
}

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Explicit transport 1D solver
///////////////////////////////////////////////////////////////////////////////

namespace {

template<Limiter L>
class TestTransport1WorkerExplicit : public ATestTransport1Worker {
    constexpr static double SLOPE_RATIO_EPS = 1e-12;

public:
    TestTransport1WorkerExplicit(size_t n_cells, double tau) : ATestTransport1Worker(n_cells, tau) {}

private:
    void impl_step() override {
        std::vector<double> fp(grid_.n_points(), 0.0);
        for (size_t i = 0; i < fp.size(); ++i) {
            // slope ratio
            double den = u_[i + 1] - u_[i];
            if (std::abs(den) < SLOPE_RATIO_EPS) {
                den = SLOPE_RATIO_EPS;
            }
            double r = (i == 0) ? 0 : (u_[i] - u_[i - 1]) / den;
            // Limiter
            double P = limiter<L>(r);

            // Low and High fluxes
            double fl = 1.0 * u_[i];
            double fh = 1.0 * (u_[i] + u_[i + 1]) / 2.0;

            // Tvd flux
            fp[i] = fl + P * (fh - fl);
        }

        u_[0] = exact_solution(grid_.point(0).x);
        u_.back() = exact_solution(grid_.point(grid_.n_points() - 1).x);
        for (size_t i = 1; i < grid_.n_points() - 1; ++i) {
            u_[i] = u_[i] - tau_ / h_ * (fp[i] - fp[i - 1]);
        }
    }
};

template<Limiter L>
class TestTransport1WorkerTheta : public ATestTransport1Worker {
    constexpr static double ITERATION_EPS = 1e-6;
    constexpr static size_t SOLVER_MAX_ITER = 100;
    constexpr static double SOLVER_EPS = 1e-8;
    constexpr static size_t MAX_ITER = 100;
    constexpr static double SLOPE_RATIO_EPS = 1e-12;

public:
    TestTransport1WorkerTheta(size_t n_cells, double tau, double theta)
        : ATestTransport1Worker(n_cells, tau), theta_(theta), E_(init_build_e()), L_(init_build_l()),
          B_(init_build_b()) {}

private:
    const double theta_;
    const LodMatrix E_;
    const LodMatrix L_;
    const LodMatrix B_;

    void impl_step() override {
        // right side
        std::vector<double> rhs = build_rhs(u_);

        double r;
        for (size_t iter = 0; iter < MAX_ITER; ++iter) {
            // left side
            LodMatrix lhs = build_lhs(u_);
            // residual and exit condition
            std::vector<double> residual = compute_residual(lhs, rhs, u_);
            r = compute_norm2(residual);
            if (r < ITERATION_EPS) {
                return;
            }
            // current iteration solution
            AmgcMatrixSolver::solve_slae(lhs.to_csr(), rhs, u_, SOLVER_MAX_ITER, SOLVER_EPS);
        }
        throw std::runtime_error("Max iteration reached with r = " + std::to_string(r));
    }

    LodMatrix init_build_e() const {
        LodMatrix ret(grid_.n_points());
        for (size_t i = 0; i < ret.n_rows(); ++i) {
            ret.set_value(i, i, 1.0);
        }
        return ret;
    };

    LodMatrix init_build_l() const {
        LodMatrix ret(grid_.n_points());
        for (size_t i = 1; i < ret.n_rows(); ++i) {
            ret.set_value(i, i, -1.0 / h_);
            ret.set_value(i, i - 1, 1.0 / h_);
        }
        return ret;
    };

    LodMatrix init_build_b() const {
        auto ret = LodMatrix::sum(1.0, -tau_ * theta_, E_, L_);
        // bc
        ret.set_unit_row(0);
        ret.set_unit_row(grid_.n_points() - 1);
        return ret;
    };

    LodMatrix build_antidiffusion(const std::vector<double>& u) const {
        LodMatrix ret(grid_.n_points());

        // F(r)
        std::vector<double> lim_values(grid_.n_points(), 0.0);
        for (size_t i = 1; i < grid_.n_points() - 1; ++i) {
            double denum = u[i + 1] - u[i];
            if (std::abs(denum) < SLOPE_RATIO_EPS) {
                denum = SLOPE_RATIO_EPS;
            }
            lim_values[i] = limiter<L>((u[i] - u[i - 1]) / denum);
        }

        for (size_t i = 1; i < grid_.n_points() - 1; ++i) {
            ret.set_value(i, i, (lim_values[i] + lim_values[i - 1]) / 2.0 / h_);
            ret.set_value(i, i - 1, -(lim_values[i - 1]) / 2.0 / h_);
            ret.set_value(i, i + 1, -(lim_values[i]) / 2.0 / h_);
        }

        return ret;
    }

    std::vector<double> build_rhs(const std::vector<double>& u_old) const {
        LodMatrix r1 = LodMatrix::sum(1.0, tau_ * (1 - theta_), E_, L_);
        LodMatrix fa = build_antidiffusion(u_old);
        LodMatrix r2 = LodMatrix::sum(1.0, tau_ * (1 - theta_), r1, fa);
        std::vector<double> ret = r2.mult_vec(u_old);
        // bc
        ret.front() = this->exact_solution(grid_.point(0).x);
        ret.back() = this->exact_solution(grid_.point(grid_.n_points() - 1).x);

        return ret;
    };

    LodMatrix build_lhs(const std::vector<double>& u) const {
        LodMatrix fa = build_antidiffusion(u);
        return LodMatrix::sum(1, -tau_ * theta_, B_, fa);
    }

    std::vector<double> compute_residual(const IMatrix& lhs, const std::vector<double>& rhs,
                                         const std::vector<double>& u) const {
        std::vector<double> ret(u.size());
        std::ranges::transform(rhs, lhs.mult_vec(u), ret.begin(), [](double f, double au) { return f - au; });
        return ret;
    }
};

} // namespace

TEST_CASE("Transport 1D solver, explicit scheme", "[transport1-fdm-theta]") {
    std::cout << std::endl << "--- cfd_test [transport1-fdm-theta] --- " << std::endl;
    // parameters
    const double tend = 0.5;
    const double tau = 1e-2;
    const double save_tau = 0.05;
    const size_t n_cells = 20;

    auto compute = [&](ATestTransport1Worker& worker, std::string out_file) {
        // non-stationary saver
        std::unique_ptr<VtkUtils::TimeSeriesWriter> writer;
        if (!out_file.empty()) {
            writer = std::make_unique<VtkUtils::TimeSeriesWriter>(out_file);
            writer->set_time_step(save_tau);
            std::string out_filename = writer->add(0);
            worker.save_vtk(out_filename);
        }

        while (worker.current_time() < tend - 1e-6) {
            // solve problem
            worker.step();

            // export solution to vtk
            if (writer) {
                if (auto fn = writer->add(worker.current_time()); !fn.empty()) {
                    worker.save_vtk(fn);
                }
            }
        };
    };

    // =========== EXPLICIT worker
    {
        // UPWIND
        TestTransport1WorkerExplicit<Limiter::UPWIND> worker(n_cells, tau);
        compute(worker, "transport1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.193719381).margin(1e-5));
    }
    {
        // MINMOD
        TestTransport1WorkerExplicit<Limiter::MINMOD> worker(n_cells, tau);
        compute(worker, "transport1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.1172775615).margin(1e-5));
    }

    // =========== THETA worker
    {
        // EXPLICIT + UPWIND
        TestTransport1WorkerTheta<Limiter::UPWIND> worker(n_cells, tau, 0.0);
        compute(worker, "transport1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.193719381).margin(1e-5));
    }

    {
        // EXPLICIT + MINMOD
        TestTransport1WorkerTheta<Limiter::MINMOD> worker(n_cells, tau, 0.0);
        compute(worker, "transport1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.1172775615).margin(1e-5));
    }

    {
        // IMPLICIT + UPWIND
        TestTransport1WorkerTheta<Limiter::UPWIND> worker(n_cells, tau, 1.0);
        compute(worker, "transport1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.2235033611).margin(1e-5));
    }

    {
        // CRANK-NICOLSON + UPWIND
        TestTransport1WorkerTheta<Limiter::UPWIND> worker(n_cells, tau, 0.5);
        compute(worker, "transport1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.2101423661).margin(1e-5));
    }

    {
        // CRANK-NICOLSON + MINMOD
        double tau1 = 1e-3;
        size_t n_cells1 = 100;
        TestTransport1WorkerTheta<Limiter::MINMOD> worker(n_cells1, tau1, 0.5);
        compute(worker, "transport1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.0235738243).margin(1e-5));
    }
}
