#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"
#include <ranges>

using namespace cfd;

namespace {

double buckley_leverett(double s) {
    s = std::min(1.0, std::max(0.0, s));
    return s * s / (s * s + (1.0 - s) * (1.0 - s));
}

double riemann_solver(double s_left, [[maybe_unused]] double s_right) {
    return s_left;
}

class ABlWorker {
public:
    virtual ~ABlWorker() {}

    static double exact_solution(Point p, [[maybe_unused]] double t) {
        return p.x < 1e-12 ? 1.0 : 0.0;
    }

    ABlWorker(std::shared_ptr<Grid1D> grid, double tau)
        : grid_(grid), tau_(tau), s_(grid->n_cells() + 2, 0), face_velocity_(grid_->n_faces(), 1.0) {}

    void step() {
        time_ += tau_;
        impl_step();
    }

    void save_vtk(const std::string& filename) {
        // save grid
        grid_->save_vtk(filename);

        // save numerical solution
        VtkUtils::add_cell_data(s_, "numerical", filename, grid_->n_cells());

        // save exact solution
        std::vector<double> exact(grid_->n_cells());
        for (size_t i = 0; i < grid_->n_cells(); ++i) {
            exact[i] = exact_solution(grid_->cell_center(i).x, time_);
        }
        VtkUtils::add_cell_data(exact, "exact", filename);
    }

    double current_time() const {
        return time_;
    }

    double compute_norm2() const {
        double norm2 = 0;
        double full_area = 0;
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double diff = s_[icell] - exact_solution(grid_->cell_center(icell), time_);
            norm2 += grid_->cell_volume(icell) * diff * diff;
            full_area += grid_->cell_volume(icell);
        }
        return std::sqrt(norm2 / full_area);
    }

    // protected:
    std::shared_ptr<IGrid> grid_;
    const double tau_;
    std::vector<double> s_;
    std::vector<double> face_velocity_;
    double time_ = 0;

private:
    virtual void impl_step() = 0;
};

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Explicit transport 1D solver
///////////////////////////////////////////////////////////////////////////////

namespace {

class BlExplicit : public ABlWorker {
public:
    BlExplicit(std::shared_ptr<Grid1D> grid, double tau) : ABlWorker(grid, tau) {}

private:
    void impl_step() override {
        // face fluxes
        std::vector<double> face_flux(grid_->n_faces());
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [cell_left, cell_right] = grid_->tab_face_cell(iface);
            double s_left = cell_left == INVALID_INDEX ? 1 : s_[cell_left];
            double s_right = cell_right == INVALID_INDEX ? 0.0 : s_[cell_right];
            // riemann solver
            double vn = face_velocity_[iface];
            double s_face = riemann_solver(s_left, s_right);
            face_flux[iface] = vn * buckley_leverett(s_face);
        }

        // explicit scheme
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [cell_left, cell_right] = grid_->tab_face_cell(iface);
            double val = tau_ * grid_->face_area(iface) * face_flux[iface];
            if (cell_left != INVALID_INDEX) {
                s_[cell_left] -= val / grid_->cell_volume(cell_left);
            }
            if (cell_right != INVALID_INDEX) {
                s_[cell_right] += val / grid_->cell_volume(cell_right);
            }
        }
    }
};

} // namespace

TEST_CASE("Transport 1D solver, explicit scheme", "[bl-explicit]") {
    std::cout << std::endl << "--- cfd_test [bl-explicit] --- " << std::endl;
    // parameters
    const double tend = 0.5;
    const double tau = 2e-2;
    const double save_tau = 0.05;
    const size_t n_cells = 20;

    auto compute = [&](ABlWorker& worker, std::string out_file) {
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
        auto grid = std::make_shared<Grid1D>(0, 1, n_cells);
        BlExplicit worker(grid, tau);
        compute(worker, "buckley_leverett1");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.6282177676).margin(1e-5));
    }
}
