#include "cfd/debug/printer.hpp"
#include "cfd/debug/tictoc.hpp"
#include "cfd/grid/regular_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"
#include "utils/vecmat.hpp"

using namespace cfd;

struct CavitySimpleWorker {
    CavitySimpleWorker(double Re, size_t n_cells, double alpha_u, double alpha_p);
    void initialize_saver(bool save_exact_fields, std::string stem, size_t iter_step);
    double set_uvp(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& p);
    double step();
    void save_current_fields(size_t iter, bool force = false);

    size_t u_size() const {
        return yf_grid_.n_points();
    }
    size_t v_size() const {
        return xf_grid_.n_points();
    }
    size_t p_size() const {
        return cc_grid_.n_points();
    }
    const std::vector<double>& pressure() const {
        return p_;
    }

    Vector top_velocity() const {
        return {1, 0};
    }

private:
    const RegularGrid2D grid_;
    const RegularGrid2D cc_grid_;
    const RegularGrid2D xf_grid_;
    const RegularGrid2D yf_grid_;
    const double hx_;
    const double hy_;
    const double Re_;
    const double alpha_u_;
    const double alpha_p_;

    std::vector<double> p_;
    std::vector<double> u_;
    std::vector<double> v_;

    const double diff_x_;
    const double diff_y_;
    const double diff_diag_;
    const double du_;

    AmgcMatrixSolver p_prime_solver_;

    CsrMatrix mat_u_;
    CsrMatrix mat_v_;
    std::vector<double> rhs_u_;
    std::vector<double> rhs_v_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_u_;
    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_v_;
    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_p_;
    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_all_;

    void assemble_p_prime_solver();
    void assemble_u_slae();
    void assemble_v_slae();

    std::vector<double> compute_u_star();
    std::vector<double> compute_v_star();
    std::vector<double> compute_p_prime(const std::vector<double>& u_star, const std::vector<double>& v_star);
    std::vector<double> compute_u_prime(const std::vector<double>& p_prime);
    std::vector<double> compute_v_prime(const std::vector<double>& p_prime);
    double u_i_j(size_t i, size_t j) const;
    double u_ip_jp(size_t i, size_t j) const;
    double v_i_j(size_t i, size_t j) const;
    double v_ip_jp(size_t i, size_t j) const;
    double p_ip_jp(size_t i, size_t j) const;
    std::vector<Vector> build_main_grid_velocity() const;
};

CavitySimpleWorker::CavitySimpleWorker(double Re, size_t n_cells, double alpha_u, double alpha_p)
    : grid_(0, 1, 0, 1, n_cells, n_cells), cc_grid_(grid_.cell_centered_grid()), xf_grid_(grid_.xface_centered_grid()),
      yf_grid_(grid_.yface_centered_grid()), hx_(1.0 / static_cast<double>(n_cells)),
      hy_(1.0 / static_cast<double>(n_cells)), Re_(Re), alpha_u_(alpha_u), alpha_p_(alpha_p),
      diff_x_(1.0 / Re_ / hx_ / hx_), diff_y_(1.0 / Re_ / hy_ / hy_), diff_diag_(2 * (diff_x_ + diff_y_)),
      du_(alpha_u_ / diff_diag_) {
    assemble_p_prime_solver();
}

void CavitySimpleWorker::initialize_saver(bool save_exact_fields, std::string stem, size_t iter_step) {
    writer_all_.reset(new VtkUtils::TimeSeriesWriter(stem));
    writer_all_->set_iter_step(iter_step);
    if (save_exact_fields) {
        writer_u_.reset(new VtkUtils::TimeSeriesWriter(stem + "-u"));
        writer_v_.reset(new VtkUtils::TimeSeriesWriter(stem + "-v"));
        writer_p_.reset(new VtkUtils::TimeSeriesWriter(stem + "-p"));

        writer_u_->set_iter_step(iter_step);
        writer_v_->set_iter_step(iter_step);
        writer_p_->set_iter_step(iter_step);
    }
};

double CavitySimpleWorker::set_uvp(const std::vector<double>& u, const std::vector<double>& v,
                                   const std::vector<double>& p) {
    u_ = u;
    v_ = v;
    p_ = p;
    assemble_u_slae();
    assemble_v_slae();

    // residuals
    auto r_u = compute_residual_vec(mat_u_, rhs_u_, u_);
    auto r_v = compute_residual_vec(mat_v_, rhs_v_, v_);
    double nrm_u = (*std::max_element(r_u.begin(), r_u.end()));
    double nrm_v = (*std::max_element(r_v.begin(), r_v.end()));

    return std::max(nrm_u, nrm_v);
};

double CavitySimpleWorker::step() {
    // Predictor step: U-star
    std::vector<double> u_star = compute_u_star();
    std::vector<double> v_star = compute_v_star();
    // Pressure correction
    std::vector<double> p_prime = compute_p_prime(u_star, v_star);
    // Velocity correction
    std::vector<double> u_prime = compute_u_prime(p_prime);
    std::vector<double> v_prime = compute_v_prime(p_prime);
    // Set final values
    std::vector<double> u_new = vector_sum(u_star, 1.0, u_prime);
    std::vector<double> v_new = vector_sum(v_star, 1.0, v_prime);
    std::vector<double> p_new = vector_sum(p_, alpha_p_, p_prime);

    return set_uvp(u_new, v_new, p_new);
}

void CavitySimpleWorker::save_current_fields(size_t iter, bool force) {
    // all data on the main grid
    if (writer_all_) {
        if (std::string filepath = writer_all_->add_iter(iter, force); !filepath.empty()) {
            grid_.save_vtk(filepath);
            VtkUtils::add_cell_data(p_, "pressure", filepath);
            VtkUtils::add_point_vector(build_main_grid_velocity(), "velocity", filepath);
        }
    }
    // pressure
    if (writer_p_) {
        if (std::string filepath = writer_p_->add_iter(iter, force); !filepath.empty()) {
            cc_grid_.save_vtk(filepath);
            VtkUtils::add_point_data(p_, "pressure", filepath);
        }
    }
    // u
    if (writer_u_) {
        if (std::string filepath = writer_u_->add_iter(iter, force); !filepath.empty()) {
            yf_grid_.save_vtk(filepath);
            VtkUtils::add_point_data(u_, "velocity-x", filepath);
        }
    }
    // v
    if (writer_v_) {
        if (std::string filepath = writer_v_->add_iter(iter, force); !filepath.empty()) {
            xf_grid_.save_vtk(filepath);
            VtkUtils::add_point_data(v_, "velocity-y", filepath);
        }
    }
}

void CavitySimpleWorker::assemble_p_prime_solver() {
    LodMatrix mat(p_size());
    for (size_t j = 0; j < cc_grid_.ny() + 1; ++j) {
        for (size_t i = 0; i < cc_grid_.nx() + 1; ++i) {
            bool is_left = (i == 0);
            bool is_right = (i == cc_grid_.nx());
            bool is_bottom = (j == 0);
            bool is_top = (j == cc_grid_.ny());

            size_t ind0 = grid_.cell_centered_grid_index_ip_jp(i, j);
            double coef_x = du_ / hx_ / hx_;
            double coef_y = du_ / hy_ / hy_;
            // x
            if (!is_right) {
                size_t ind1 = grid_.cell_centered_grid_index_ip_jp(i + 1, j);
                mat.add_value(ind0, ind0, coef_x);
                mat.add_value(ind0, ind1, -coef_x);
            }
            if (!is_left) {
                size_t ind1 = grid_.cell_centered_grid_index_ip_jp(i - 1, j);
                mat.add_value(ind0, ind0, coef_x);
                mat.add_value(ind0, ind1, -coef_x);
            }
            // y
            if (!is_top) {
                size_t ind1 = grid_.cell_centered_grid_index_ip_jp(i, j + 1);
                mat.add_value(ind0, ind0, coef_y);
                mat.add_value(ind0, ind1, -coef_y);
            }
            if (!is_bottom) {
                size_t ind1 = grid_.cell_centered_grid_index_ip_jp(i, j - 1);
                mat.add_value(ind0, ind0, coef_y);
                mat.add_value(ind0, ind1, -coef_y);
            }
        }
    }
    // p[0] = 0
    mat.set_unit_row(0);

    p_prime_solver_.set_matrix(mat.to_csr());
}

void CavitySimpleWorker::assemble_u_slae() {
    rhs_u_.resize(u_.size());
    std::fill(rhs_u_.begin(), rhs_u_.end(), 0.0);
    LodMatrix mat(u_.size());

    auto add_to_mat = [&](size_t row_index, std::array<size_t, 2> ij_col, double value) {
        if (ij_col[1] == grid_.ny()) {
            // ghost index => top boundary condition: u = u_top
            size_t ind1 = grid_.yface_grid_index_i_jp(ij_col[0], ij_col[1] - 1);
            mat.add_value(row_index, ind1, -value);
            rhs_u_[row_index] -= 2.0 * value * top_velocity().x;
        } else if (ij_col[1] == (size_t)-1) {
            // ghost index => bottom boundary condition: u = 0
            size_t ind1 = grid_.yface_grid_index_i_jp(ij_col[0], ij_col[1] + 1);
            mat.add_value(row_index, ind1, -value);
        } else {
            size_t ind1 = grid_.yface_grid_index_i_jp(ij_col[0], ij_col[1]);
            mat.add_value(row_index, ind1, value);
        }
    };

    // left/right boundary: u = 0
    for (size_t j = 0; j < grid_.ny(); ++j) {
        size_t index_left = grid_.yface_grid_index_i_jp(0, j);
        add_to_mat(index_left, {0, j}, 1.0);
        rhs_u_[index_left] = 0.0;

        size_t index_right = grid_.yface_grid_index_i_jp(grid_.nx(), j);
        add_to_mat(index_right, {grid_.nx(), j}, 1.0);
        rhs_u_[index_right] = 0.0;
    }

    // internal
    for (size_t j = 0; j < grid_.ny(); ++j) {
        for (size_t i = 1; i < grid_.nx(); ++i) {
            size_t row_index = grid_.yface_grid_index_i_jp(i, j); //[i, j+1/2]

            double u0_plus = u_ip_jp(i, j);      //_u[i+1/2, j+1/2]
            double u0_minus = u_ip_jp(i - 1, j); //_u[i-1/2, j+1/2]
            double v0_plus = v_i_j(i, j + 1);    //_v[i,j+1]
            double v0_minus = v_i_j(i, j);       //_v[i,j]

            // diagonal diffusion
            add_to_mat(row_index, {i, j}, diff_diag_ / alpha_u_);
            // offdiagonal diffusion
            add_to_mat(row_index, {i + 1, j}, -diff_x_);
            add_to_mat(row_index, {i - 1, j}, -diff_x_);
            add_to_mat(row_index, {i, j + 1}, -diff_y_);
            add_to_mat(row_index, {i, j - 1}, -diff_y_);
            // + d(u0*u)/dx
            add_to_mat(row_index, {i + 1, j}, 1.0 / 2.0 / hx_ * u0_plus);
            add_to_mat(row_index, {i - 1, j}, -1.0 / 2.0 / hx_ * u0_minus);
            // + d(v0*u)/dy
            add_to_mat(row_index, {i, j + 1}, 1.0 / 2.0 / hy_ * v0_plus);
            add_to_mat(row_index, {i, j - 1}, -1.0 / 2.0 / hy_ * v0_minus);
            // = diagonal diffusion relaxation
            rhs_u_[row_index] += (1.0 - alpha_u_) / alpha_u_ * diff_diag_ * u_[row_index];
            // - dp/dx
            rhs_u_[row_index] -= 1.0 / hx_ * (p_ip_jp(i, j) - p_ip_jp(i - 1, j));
        }
    }
    mat_u_ = mat.to_csr();
}

void CavitySimpleWorker::assemble_v_slae() {
    rhs_v_.resize(v_.size());
    std::fill(rhs_v_.begin(), rhs_v_.end(), 0.0);
    LodMatrix mat(v_.size());

    auto add_to_mat = [&](size_t row_index, std::array<size_t, 2> ij_col, double value) {
        if (ij_col[0] == (size_t)-1) {
            // left boundary condition: v = 0
            size_t ind1 = grid_.xface_grid_index_ip_j(ij_col[0] + 1, ij_col[1]);
            mat.add_value(row_index, ind1, -value);
        } else if (ij_col[0] == grid_.nx()) {
            // right boundary condition: v = 0
            size_t ind1 = grid_.xface_grid_index_ip_j(ij_col[0] - 1, ij_col[1]);
            mat.add_value(row_index, ind1, -value);
        } else {
            // internal
            size_t ind1 = grid_.xface_grid_index_ip_j(ij_col[0], ij_col[1]);
            mat.add_value(row_index, ind1, value);
        }
    };

    // === top/bottom boundaries
    for (size_t i = 0; i < grid_.nx(); ++i) {
        // top boundary condition: v = 0;
        size_t index_top = grid_.xface_grid_index_ip_j(i, grid_.ny());
        add_to_mat(index_top, {i, grid_.ny()}, 1.0);
        rhs_v_[index_top] = 0.0;
        // bottom boundary condition: v = 0;
        size_t index_bottom = grid_.xface_grid_index_ip_j(i, 0);
        add_to_mat(index_bottom, {i, 0}, 1.0);
        rhs_v_[index_bottom] = 0.0;
    }

    // === internal
    for (size_t j = 1; j < grid_.ny(); ++j) {
        for (size_t i = 0; i < grid_.nx(); ++i) {
            size_t row_index = grid_.xface_grid_index_ip_j(i, j); //[i+1/2, j] linear index in v grid

            double u0_plus = u_i_j(i + 1, j);    // _u[i+1, j]
            double u0_minus = u_i_j(i, j);       // _u[i, j]
            double v0_plus = v_ip_jp(i, j);      // _v[i+1/2, j+1/2]
            double v0_minus = v_ip_jp(i, j - 1); // _v[i+1/2, j-1/2]

            // diagonal diffusion
            add_to_mat(row_index, {i, j}, diff_diag_ / alpha_u_);
            // offdiagonal diffusion
            add_to_mat(row_index, {i + 1, j}, -diff_x_);
            add_to_mat(row_index, {i - 1, j}, -diff_x_);
            add_to_mat(row_index, {i, j + 1}, -diff_y_);
            add_to_mat(row_index, {i, j - 1}, -diff_y_);
            // + d (u0*v) / dx
            add_to_mat(row_index, {i + 1, j}, 1.0 / 2.0 / hx_ * u0_plus);
            add_to_mat(row_index, {i - 1, j}, -1.0 / 2.0 / hx_ * u0_minus);
            // + d (v0*v) / dy
            add_to_mat(row_index, {i, j + 1}, 1.0 / 2.0 / hy_ * v0_plus);
            add_to_mat(row_index, {i, j - 1}, -1.0 / 2.0 / hy_ * v0_minus);
            // = diagonal diffusion relaxation
            rhs_v_[row_index] += (1.0 - alpha_u_) / alpha_u_ * diff_diag_ * v_[row_index];
            // - dp/dy
            rhs_v_[row_index] -= 1.0 * (p_ip_jp(i, j) - p_ip_jp(i, j - 1)) / hy_;
        }
    }
    mat_v_ = mat.to_csr();
}

std::vector<double> CavitySimpleWorker::compute_u_star() {
    std::vector<double> u_star(u_);
    AmgcMatrixSolver::solve_slae(mat_u_, rhs_u_, u_star);
    return u_star;
}

std::vector<double> CavitySimpleWorker::compute_v_star() {
    std::vector<double> v_star(v_);
    AmgcMatrixSolver::solve_slae(mat_v_, rhs_v_, v_star);
    return v_star;
}

std::vector<double> CavitySimpleWorker::compute_p_prime(const std::vector<double>& u_star,
                                                        const std::vector<double>& v_star) {
    // assemble rhs = du/dx + dv/dy
    std::vector<double> rhs(p_.size(), 0.0);
    for (size_t i = 0; i < grid_.nx(); ++i) {
        for (size_t j = 0; j < grid_.ny(); ++j) {
            size_t ind0 = grid_.cell_centered_grid_index_ip_jp(i, j);
            size_t ind_left = grid_.yface_grid_index_i_jp(i, j);
            size_t ind_right = grid_.yface_grid_index_i_jp(i + 1, j);
            size_t ind_bot = grid_.xface_grid_index_ip_j(i, j);
            size_t ind_top = grid_.xface_grid_index_ip_j(i, j + 1);
            rhs[ind0] = -(u_star[ind_right] - u_star[ind_left]) / hx_ - (v_star[ind_top] - v_star[ind_bot]) / hy_;
        }
    }
    // p[0] = 0
    rhs[0] = 0;

    std::vector<double> p_prime;
    p_prime_solver_.solve(rhs, p_prime);
    return p_prime;
}

std::vector<double> CavitySimpleWorker::compute_u_prime(const std::vector<double>& p_prime) {
    std::vector<double> u_prime(u_.size(), 0.0);
    for (size_t i = 1; i < grid_.nx(); ++i) {
        for (size_t j = 0; j < grid_.ny(); ++j) {
            size_t ind0 = grid_.yface_grid_index_i_jp(i, j);
            size_t ind_plus = grid_.cell_centered_grid_index_ip_jp(i, j);
            size_t ind_minus = grid_.cell_centered_grid_index_ip_jp(i - 1, j);
            u_prime[ind0] = -du_ * (p_prime[ind_plus] - p_prime[ind_minus]) / hx_;
        }
    }
    return u_prime;
}

std::vector<double> CavitySimpleWorker::compute_v_prime(const std::vector<double>& p_prime) {
    std::vector<double> v_prime(v_.size(), 0.0);
    for (size_t i = 0; i < grid_.nx(); ++i) {
        for (size_t j = 1; j < grid_.ny(); ++j) {
            size_t ind0 = grid_.xface_grid_index_ip_j(i, j);
            size_t ind_plus = grid_.cell_centered_grid_index_ip_jp(i, j);
            size_t ind_minus = grid_.cell_centered_grid_index_ip_jp(i, j - 1);
            v_prime[ind0] = -du_ * (p_prime[ind_plus] - p_prime[ind_minus]) / hy_;
        }
    }
    return v_prime;
}

double CavitySimpleWorker::u_i_j(size_t i, size_t j) const {
    size_t ind0 = grid_.yface_grid_index_i_jp(i, j);
    size_t ind1 = grid_.yface_grid_index_i_jp(i, j - 1);
    return (u_[ind0] + u_[ind1]) / 2.0;
}
double CavitySimpleWorker::u_ip_jp(size_t i, size_t j) const {
    size_t ind0 = grid_.yface_grid_index_i_jp(i, j);
    size_t ind1 = grid_.yface_grid_index_i_jp(i + 1, j);
    return (u_[ind0] + u_[ind1]) / 2.0;
}

double CavitySimpleWorker::v_i_j(size_t i, size_t j) const {
    size_t ind0 = grid_.xface_grid_index_ip_j(i, j);
    size_t ind1 = grid_.xface_grid_index_ip_j(i - 1, j);
    return (v_[ind0] + v_[ind1]) / 2.0;
}
double CavitySimpleWorker::v_ip_jp(size_t i, size_t j) const {
    size_t ind0 = grid_.xface_grid_index_ip_j(i, j);
    size_t ind1 = grid_.xface_grid_index_ip_j(i, j + 1);
    return (v_[ind0] + v_[ind1]) / 2.0;
}

double CavitySimpleWorker::p_ip_jp(size_t i, size_t j) const {
    return p_[grid_.cell_centered_grid_index_ip_jp(i, j)];
}

std::vector<Vector> CavitySimpleWorker::build_main_grid_velocity() const {
    std::vector<Vector> ret(grid_.n_points());
    // boundary
    for (size_t j = 0; j < grid_.ny() + 1; ++j) {
        // left boundary
        size_t ind_left = grid_.to_linear_point_index({0, j});
        ret[ind_left] = Vector(0, 0);
        // right boundary
        size_t ind_right = grid_.to_linear_point_index({grid_.nx(), j});
        ret[ind_right] = Vector(0, 0);
    }
    for (size_t i = 0; i < grid_.nx() + 1; ++i) {
        // bottom boundary
        size_t ind_bot = grid_.to_linear_point_index({i, 0});
        ret[ind_bot] = 0;
        // top boundary
        size_t ind_top = grid_.to_linear_point_index({i, grid_.ny()});
        ret[ind_top] = top_velocity();
    }

    // internal
    for (size_t j = 1; j < grid_.ny(); ++j) {
        for (size_t i = 1; i < grid_.nx(); ++i) {
            size_t ind = grid_.to_linear_point_index({i, j});
            size_t ind_top = grid_.yface_grid_index_i_jp(i, j);
            size_t ind_bot = grid_.yface_grid_index_i_jp(i, j - 1);
            size_t ind_left = grid_.xface_grid_index_ip_j(i - 1, j);
            size_t ind_right = grid_.xface_grid_index_ip_j(i, j);
            ret[ind] = Vector(0.5 * (u_[ind_top] + u_[ind_bot]), 0.5 * (v_[ind_left] + v_[ind_right]));
        }
    }
    return ret;
}

TEST_CASE("Cavity, SIMPLE fdm algorithm", "[cavity-fdm-simple]") {
    std::cout << std::endl << "--- cfd_test [cavity-fdm-simple] --- " << std::endl;

    // problem parameters
    double Re = 100;
    double alpha_u = 0.8;
    double alpha_p = 0.3;
    size_t n_cells = 30;
    size_t max_it = 1000;
    double eps = 1e-0;

    // worker initialization
    CavitySimpleWorker worker(Re, n_cells, alpha_u, alpha_p);
    worker.initialize_saver(false, "cavity-fdm", 5);

    // initial condition
    std::vector<double> u_init(worker.u_size(), 0.0);
    std::vector<double> v_init(worker.v_size(), 0.0);
    std::vector<double> p_init(worker.p_size(), 0.0);
    worker.set_uvp(u_init, v_init, p_init);
    worker.save_current_fields(0);

    // iterations loop
    size_t it;
    for (it = 1; it < max_it; ++it) {
        double nrm = worker.step();

        // print norm and pressure value at the top-right corner
        std::cout << it << " " << nrm << " " << worker.pressure().back() << std::endl;

        // export solution to vtk
        worker.save_current_fields(it);

        // break if residual is low enough
        if (nrm < eps) {
            break;
        }
    }
    worker.save_current_fields(it, true);

    CHECK(it == 9);
}
