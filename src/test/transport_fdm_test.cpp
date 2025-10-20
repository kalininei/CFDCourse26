#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"

using namespace cfd;

class ATestTransport1Worker {
public:
    virtual ~ATestTransport1Worker() {}

    static double init_solution(double x) {
        // return (x >= -0.1 && x <= 0.1) ? 1.0 : 0.0;
        constexpr double sigma = 0.1;
        return exp(-x * x / sigma / sigma);
    }
    double exact_solution(double x) {
        return init_solution(x - time_);
    }

    ATestTransport1Worker(size_t n_cells)
        : grid_(0, 1, n_cells), h_(1.0 / static_cast<double>(n_cells)), u_(grid_.n_points()) {
        for (size_t i = 0; i < grid_.n_points(); ++i) {
            u_[i] = init_solution(grid_.point(i).x);
        }
    }

    double step(double tau) {
        time_ += tau;
        impl_step(tau);
        return compute_norm2();
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
        // weights
        std::vector<double> w(grid_.n_points(), h_);
        w[0] = w[grid_.n_points() - 1] = h_ / 2;

        // sum
        double sum = 0;
        for (size_t i = 0; i < grid_.n_points(); ++i) {
            double diff = u_[i] - exact_solution(grid_.point(i).x);
            sum += w[i] * diff * diff;
        }

        double len = grid_.point(grid_.n_points() - 1).x - grid_.point(0).x;
        return std::sqrt(sum / len);
    }

protected:
    Grid1D grid_;
    double h_;
    std::vector<double> u_;
    double time_ = 0;

private:
    virtual void impl_step(double tau) = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Explicit transport 1D solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport1WorkerExplicit : public ATestTransport1Worker {
public:
    TestTransport1WorkerExplicit(size_t n_cells) : ATestTransport1Worker(n_cells) {}

private:
    void impl_step(double tau) override {
        std::vector<double> u_old(u_);
        u_[0] = exact_solution(grid_.point(0).x);
        for (size_t i = 1; i < grid_.n_points(); ++i) {
            u_[i] = u_old[i] - tau / h_ * (u_old[i] - u_old[i - 1]);
        }
    }
};

TEST_CASE("Transport 1D solver, explicit", "[transport1-fdm-explicit]") {
    std::cout << std::endl << "--- cfd_test [transport1-fdm-explicit] --- " << std::endl;
    // parameters
    const double tend = 0.5;
    const double V = 1.0;
    const double L = 1.0;
    size_t n_cells = 100;
    double Cu = 0.9;
    double h = L / static_cast<double>(n_cells);
    double tau = Cu * h / V;

    for (size_t n_cells2 : std::vector<size_t>{10, 50, 100, 300, 500, 1000, 3000, 5000, 10000, 30000, 50000}) {
        n_cells = n_cells2;
        tau = 1e-5;
        // solver
        TestTransport1WorkerExplicit worker(n_cells);

        // saver
        // VtkUtils::TimeSeriesWriter writer("transport1-explicit");
        // std::string out_filename = writer.add(0);
        // worker.save_vtk(out_filename);

        double norm = 0;
        while (worker.current_time() < tend - 1e-6) {
            // solve problem
            norm += tau * worker.step(tau);
            // export solution to vtk
            // out_filename = writer.add(worker.current_time());
            // worker.save_vtk(out_filename);
        };
        std::cout << n_cells << " " << norm << std::endl;
    }
    // std::cout << 1.0/tau << " " << norm << std::endl;
    // CHECK(norm == Approx(0.0138123932).margin(1e-5));
}

///////////////////////////////////////////////////////////////////////////////
// Implicit transport 1D solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport1WorkerImplicit : public ATestTransport1Worker {
public:
    TestTransport1WorkerImplicit(size_t n_cells) : ATestTransport1Worker(n_cells) {}

private:
    void impl_step(double tau) override {
        AmgcMatrixSolver& slv = build_solver(tau);
        std::vector<double> rhs = build_rhs(tau);
        slv.solve(rhs, u_);
    }

    AmgcMatrixSolver _solver;
    double _last_used_tau = 0;

    AmgcMatrixSolver& build_solver(double tau) {
        if (std::abs(_last_used_tau - tau) > 1e-12) {
            CsrMatrix mat = build_lhs(tau);
            _solver.set_matrix(mat);
            _last_used_tau = tau;
        }
        return _solver;
    }

    virtual CsrMatrix build_lhs(double tau) {
        LodMatrix mat(u_.size());
        mat.set_value(0, 0, 1.0);
        mat.set_value(u_.size() - 1, u_.size() - 1, 1.0);
        double diag = 1.0 + tau / h_;
        double nondiag = -tau / h_;
        for (size_t i = 1; i < u_.size() - 1; ++i) {
            mat.set_value(i, i, diag);
            mat.set_value(i, i - 1, nondiag);
        }
        return mat.to_csr();
    }

    virtual std::vector<double> build_rhs([[maybe_unused]] double tau) {
        std::vector<double> rhs(u_);
        rhs[0] = exact_solution(grid_.point(0).x);
        rhs.back() = exact_solution(grid_.point(grid_.n_points() - 1).x);
        return rhs;
    }
};

TEST_CASE("Transport 1D solver, implicit", "[transport1-fdm-implicit]") {
    std::cout << std::endl << "--- cfd_test [transport1-fdm-implicit] --- " << std::endl;
    const double tend = 0.5;
    const double V = 1.0;
    const double L = 1.0;
    size_t n_cells = 100;
    double Cu = 0.9;
    double h = L / static_cast<double>(n_cells);
    double tau = Cu * h / V;

    for (size_t n_cells2 : std::vector<size_t>{10, 50, 100, 300, 500, 1000, 3000, 5000, 10000, 30000, 50000}) {
        n_cells = n_cells2;
        tau = 1e-5;
        TestTransport1WorkerImplicit worker(n_cells);

        // VtkUtils::TimeSeriesWriter writer("transport1-implicit");
        // worker.save_vtk(writer.add(0));
        double norm = 0.0;
        while (worker.current_time() < tend - 1e-6) {
            norm += tau * worker.step(tau);
            // worker.save_vtk(writer.add(worker.current_time()));
        };
        std::cout << n_cells << " " << norm << std::endl;
    }
    // std::cout << 1.0/tau << " " << norm << std::endl;
    // CHECK(norm == Approx(0.137664).margin(1e-5));
}

///////////////////////////////////////////////////////////////////////////////
// Crank-Nicolson transport 1D solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport1WorkerCN : public TestTransport1WorkerImplicit {
public:
    TestTransport1WorkerCN(size_t n_cells) : TestTransport1WorkerImplicit(n_cells) {}

private:
    CsrMatrix build_lhs(double tau) override {
        LodMatrix mat(u_.size());
        mat.set_value(0, 0, 1.0);
        mat.set_value(u_.size() - 1, u_.size() - 1, 1.0);
        double diag = 1.0 + 0.5 * tau / h_;
        double nondiag = -0.5 * tau / h_;
        for (size_t i = 1; i < u_.size() - 1; ++i) {
            mat.set_value(i, i, diag);
            mat.set_value(i, i - 1, nondiag);
        }
        return mat.to_csr();
    }

    std::vector<double> build_rhs(double tau) override {
        std::vector<double> rhs(u_);
        rhs[0] = exact_solution(grid_.point(0).x);
        rhs.back() = exact_solution(grid_.point(grid_.n_points() - 1).x);
        for (size_t i = 1; i < rhs.size() - 1; ++i) {
            rhs[i] -= 0.5 * tau / h_ * (u_[i] - u_[i - 1]);
        }
        return rhs;
    }
};

TEST_CASE("Transport 1D solver, Crank-Nicolson", "[transport1-fdm-cn]") {
    std::cout << std::endl << "--- cfd_test [transport1-fdm-cn] --- " << std::endl;
    const double tend = 0.5;
    const double V = 1.0;
    const double L = 1.0;
    size_t n_cells = 100;
    double Cu = 0.9;

    double h = L / (double)n_cells;
    double tau = Cu * h / V;

    for (size_t n_cells2 : std::vector<size_t>{10, 50, 100, 300, 500, 1000, 3000, 5000, 10000, 30000, 50000}) {
        n_cells = n_cells2;
        tau = 1e-5;
        TestTransport1WorkerCN worker(n_cells);

        // VtkUtils::TimeSeriesWriter writer("transport1-cn");
        // worker.save_vtk(writer.add(0));

        double norm = 0;
        while (worker.current_time() < tend - 1e-6) {
            norm += tau * worker.step(tau);
            // worker.save_vtk(writer.add(worker.current_time()));
        };
        std::cout << n_cells << " " << norm << std::endl;
    }
    // std::cout << 1.0/tau << " " << norm << std::endl;
    // CHECK(norm == Approx(0.0937748).margin(1e-5));
}

///////////////////////////////////////////////////////////////////////////////
// Explicit tvd-transport 1D solver
///////////////////////////////////////////////////////////////////////////////

class TestTransport1WorkerTvdExplicit : public ATestTransport1Worker {
public:
    TestTransport1WorkerTvdExplicit(size_t n_cells) : ATestTransport1Worker(n_cells) {}

private:
    void impl_step(double tau) override {
        std::vector<double> fp(grid_.n_points() - 1);
        for (size_t i = 0; i < fp.size(); ++i) {
            double fl = 1.0 * u_[i];
            double fh = 1.0 * (u_[i] + u_[i + 1]) / 2.0;

            double den = u_[i + 1] - u_[i];
            if (std::abs(den) < 1e-16) {
                den += 1e-8;
            }
            double r = (i == 0) ? 1.0 : (u_[i] - u_[i - 1]) / den;

            double P = 0;
            if (r > 0) {
                P = std::max(0.0, std::min(1.0, r)); // minmod
                // P = (r < 1.0) ? std::min(1.0, 2 * r) : std::min(2.0, r);   // Superbee
            }
            fp[i] = fl + P * (fh - fl);
        }

        u_[0] = exact_solution(grid_.point(0).x);
        u_.back() = exact_solution(grid_.point(grid_.n_points() - 1).x);
        for (size_t i = 1; i < grid_.n_points() - 1; ++i) {
            u_[i] = u_[i] - tau / h_ * (fp[i] - fp[i - 1]);
        }
    }
};

TEST_CASE("Transport 1D solver, tvd-explicit", "[transport1-fdm-tvd-explicit]") {
    std::cout << std::endl << "--- cfd_test [transport1-fdm-tvd-explicit] --- " << std::endl;
    // parameters
    const double tend = 0.5;
    const double V = 1.0;
    const double L = 1.0;
    size_t n_cells = 100;
    double Cu = 0.9;
    double h = L / static_cast<double>(n_cells);
    double tau = Cu * h / V;

    for (size_t n_cells2 : std::vector<size_t>{10, 50, 100, 300, 500, 1000, 3000, 5000, 10000, 30000, 50000}) {
        n_cells = n_cells2;
        tau = 1e-5;
        // solver
        TestTransport1WorkerTvdExplicit worker(n_cells);

        // saver
        // VtkUtils::TimeSeriesWriter writer("transport1-explicit");
        // std::string out_filename = writer.add(0);
        // worker.save_vtk(out_filename);

        double norm = 0;
        while (worker.current_time() < tend - 1e-6) {
            // solve problem
            norm += tau * worker.step(tau);
            // export solution to vtk
            // out_filename = writer.add(worker.current_time());
            // worker.save_vtk(out_filename);
        };
        std::cout << n_cells << " " << norm << std::endl;
    }
    // std::cout << 1.0/tau << " " << norm << std::endl;
    // CHECK(norm == Approx(0.0138123932).margin(1e-5));
}
