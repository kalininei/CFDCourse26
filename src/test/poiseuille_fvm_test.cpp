#include "cfd/debug/printer.hpp"
#include "cfd/debug/saver.hpp"
#include "cfd/debug/tictoc.hpp"
#include "cfd/fvm/fvm_dfdn.hpp"
#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/fvm/fvm_gradient.hpp"
#include "cfd/grid/regular_grid2d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/matrix_iter.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"
#include "utils/filesystem.hpp"
#include "utils/vecmat.hpp"
#include <iomanip>
#include <list>

using namespace cfd;

namespace {

double input_profile(double y) {
    return -6.0 * (y - 0.5) * (y + 0.5);
};

struct PoiseuilleFvmWorker {
    PoiseuilleFvmWorker(const IGrid& grid, double Re, double alpha_u, double alpha_p);
    void initialize_saver(std::string stem, size_t save_iter_step);

    double step();
    void save_current_fields(size_t iter) const;
    void to_next_time_step() {};

    size_t vec_size() const {
        return collocations_.size();
    }

private:
    const IGrid& grid_;
    const double Re_;
    const double alpha_u_;
    const double alpha_p_;
    const FvmExtendedCollocations collocations_;
    const FvmFacesDn df_dn_computer_;
    const std::shared_ptr<IFvmCellGradient> grad_computer_;

    struct BoundaryInfo {
        std::vector<size_t> all;
        std::vector<size_t> wall;
        std::vector<size_t> input;
        std::vector<size_t> output;
    };
    BoundaryInfo boundary_info_;

    std::vector<double> p_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> du_;
    std::vector<double> dv_;
    std::vector<double> un_face_;
    std::vector<double> dn_face_;
    std::vector<double> dpdn_face_;
    std::vector<Vector> grad_p_;

    CsrMatrix mat_u_;
    CsrMatrix mat_v_;
    std::vector<double> rhs_u_;
    std::vector<double> rhs_v_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_;

    CsrMatrix assemble_u_matrix() const;
    void prepare_slae(const CsrMatrix& A);
    std::vector<double> compute_u_tilde(const CsrMatrix& A, const std::vector<double>& d,
                                        const std::vector<double>& u_star, const std::vector<double>& u) const;
    std::vector<double> compute_p(const std::vector<double>& un_face) const;
    std::vector<double> compute_un_face(const std::vector<double>& x, const std::vector<double>& y) const;
    std::vector<double> compute_dn_face(const std::vector<double>& x, const std::vector<double>& y) const;
    std::vector<double> compute_un_star_face_rhie_chow(const std::vector<double>& u_star,
                                                       const std::vector<double>& v_star) const;

    double to_next_iteration();
    void gather_boundary_collocations();
};

PoiseuilleFvmWorker::PoiseuilleFvmWorker(const IGrid& grid, double Re, double alpha_u, double alpha_p)
    : grid_(grid),
      Re_(Re),
      alpha_u_(alpha_u),
      alpha_p_(alpha_p),
      collocations_(grid),
      df_dn_computer_(grid, collocations_),
      grad_computer_(std::make_shared<GaussLinearFvmCellGradient>(grid, collocations_)) {

    gather_boundary_collocations();

    rhs_u_ = std::vector<double>(vec_size(), 0);
    rhs_v_ = std::vector<double>(vec_size(), 0);
    du_ = std::vector<double>(vec_size(), 0);
    dv_ = std::vector<double>(vec_size(), 0);
    u_ = std::vector<double>(vec_size(), 1);
    v_ = std::vector<double>(vec_size(), 0);
    p_ = std::vector<double>(vec_size(), 0);
    un_face_ = std::vector<double>(grid_.n_faces(), 0);

    dn_face_ = std::vector<double>(grid_.n_faces(), 0);
    dpdn_face_ = std::vector<double>(grid_.n_faces(), 0);
    grad_p_ = std::vector<Vector>(grid_.n_cells(), Vector{0, 0, 0});

    to_next_iteration();
}

void PoiseuilleFvmWorker::gather_boundary_collocations() {
    double xmin = grid_.point(0).x;
    double xmax = grid_.point(0).x;
    double ymin = grid_.point(0).y;
    double ymax = grid_.point(0).y;
    for (size_t i = 1; i < grid_.n_points(); ++i) {
        Point p = grid_.point(i);
        xmin = std::min(xmin, p.x);
        ymin = std::min(ymin, p.y);
        xmax = std::max(xmax, p.x);
        ymax = std::max(ymax, p.y);
    }

    BoundaryInfo& bi = boundary_info_;
    for (size_t icolloc: collocations_.face_collocations) {
        size_t iface = collocations_.face_index(icolloc);
        bi.all.push_back(icolloc);
        Point fc = grid_.face_center(iface);
        if (std::abs(fc.y - ymin) < 1e-6) {
            bi.wall.push_back(icolloc);
        } else if (std::abs(fc.y - ymax) < 1e-6) {
            bi.wall.push_back(icolloc);
        } else if (std::abs(fc.x - xmin) < 1e-6) {
            bi.input.push_back(icolloc);
        } else if (std::abs(fc.x - xmax) < 1e-6) {
            bi.output.push_back(icolloc);
        }
    }
}

void PoiseuilleFvmWorker::initialize_saver(std::string stem, size_t save_iter_step) {
    writer_.reset(new VtkUtils::TimeSeriesWriter(stem));
    writer_->set_iter_step(save_iter_step);
};

double PoiseuilleFvmWorker::to_next_iteration() {
    CsrMatrix A = assemble_u_matrix();
    prepare_slae(A);

    // residual vectors
    std::vector<double> res_u = compute_residual_vec(mat_u_, rhs_u_, u_);
    std::vector<double> res_v = compute_residual_vec(mat_v_, rhs_v_, v_);
    // norm
    double res = 0;
    for (size_t icell = 0; icell < grid_.n_cells(); ++icell) {
        res = std::max(res, std::max(res_u[icell], res_v[icell]) / grid_.cell_volume(icell));
    }
    return res;
};

std::vector<double> PoiseuilleFvmWorker::compute_un_star_face_rhie_chow(const std::vector<double>& u_star,
                                                                        const std::vector<double>& v_star) const {

    std::vector<double> ret(grid_.n_faces());

    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [ci, cj] = grid_.tab_face_cell(iface);
        if (ci == INVALID_INDEX || cj == INVALID_INDEX) {
            size_t icolloc = collocations_.boundary_colloc(iface);
            ret[iface] = dot_product(Vector{u_star[icolloc], v_star[icolloc]}, grid_.face_normal(iface));
        } else {
            // Rhie-Chow interpolation
            Vector normal = grid_.face_normal(iface);
            Vector uvec_i = Vector(u_star[ci], v_star[ci]);
            Vector uvec_j = Vector(u_star[cj], v_star[cj]);

            double ustar_i = dot_product(uvec_i, normal);
            double ustar_j = dot_product(uvec_j, normal);
            double d_dpdn_i = dot_product(Vector{grad_p_[ci].x * du_[ci], grad_p_[ci].y * dv_[ci]}, normal);
            double d_dpdn_j = dot_product(Vector{grad_p_[cj].x * du_[cj], grad_p_[cj].y * dv_[cj]}, normal);
            double d_dpdn_ij = dn_face_[iface] * dpdn_face_[iface];

            ret[iface] = 0.5 * (ustar_i + ustar_j) + 0.5 * (d_dpdn_i + d_dpdn_j - 2 * d_dpdn_ij);
        }
    }

    return ret;
}

double PoiseuilleFvmWorker::step() {
    // 1. Predictor step: U-star
    std::vector<double> u_star = u_, v_star = v_;
    AmgcMatrixSolver::solve_slae(mat_u_, rhs_u_, u_star);
    AmgcMatrixSolver::solve_slae(mat_v_, rhs_v_, v_star);
    std::vector<double> un_star_face = compute_un_star_face_rhie_chow(u_star, v_star);

    // 2. Pressure correction
    std::vector<double> p_prime = compute_p(un_star_face);
    std::vector<Vector> grad_p_prime = grad_computer_->compute(p_prime);
    std::vector<double> dpprime_dn_face = df_dn_computer_.compute(p_prime);
    // 3. Velocity correction
    std::vector<double> u_prime(vec_size(), 0.0), v_prime(vec_size(), 0.0), un_prime_face(grid_.n_faces());
    for (size_t i = 0; i < grid_.n_cells(); ++i) {
        u_prime[i] = -du_[i] * grad_p_prime[i].x;
        v_prime[i] = -dv_[i] * grad_p_prime[i].y;
    }
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        un_prime_face[iface] = -dn_face_[iface] * dpprime_dn_face[iface];
    }
    // 4. Set final values
    u_ = vector_sum(u_star, 1.0, u_prime);
    v_ = vector_sum(v_star, 1.0, v_prime);
    un_face_ = vector_sum(un_star_face, 1.0, un_prime_face);
    p_ = vector_sum(p_, alpha_p_, p_prime);
    grad_p_ = vector_sum(grad_p_, alpha_p_, grad_p_prime);
    dpdn_face_ = vector_sum(dpdn_face_, alpha_p_, dpprime_dn_face);

    return to_next_iteration();
}

void PoiseuilleFvmWorker::save_current_fields(size_t iter) const {
    if (writer_) {
        std::string filepath = writer_->add_iter(iter);
        if (filepath.empty()) {
            return;
        }
        grid_.save_vtk(filepath);
        VtkUtils::add_cell_data(p_, "pressure", filepath, grid_.n_cells());
        VtkUtils::add_cell_vector(u_, v_, "velocity", filepath, grid_.n_cells());
    }
}

std::vector<double> PoiseuilleFvmWorker::compute_dn_face(const std::vector<double>& x,
                                                         const std::vector<double>& y) const {
    std::vector<double> ret(grid_.n_faces(), 0.0);
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [ci, cj] = collocations_.tab_face_colloc(iface);
        Vector normal = grid_.face_normal(iface);
        Vector u_face;
        if (collocations_.is_boundary_colloc(ci)) {
            u_face = Vector{x[cj], y[cj]};
        } else if (collocations_.is_boundary_colloc(cj)) {
            u_face = Vector{x[ci], y[ci]};
        } else {
            u_face = Vector{0.5 * (x[ci] + x[cj]), 0.5 * (y[ci] + y[cj])};
        }
        // d_nn = d_xx * nx^2 + d_yy * ny^2, d_ns = 0 assuming d_xx = d_yy
        ret[iface] = u_face.x * normal.x * normal.x + u_face.y * normal.y * normal.y;
    }
    return ret;
}

CsrMatrix PoiseuilleFvmWorker::assemble_u_matrix() const {
    LodMatrix mat(vec_size());
    double coef_diff = 1.0 / Re_;

    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [negative_colloc, positive_colloc] = collocations_.tab_face_colloc(iface);
        double area = grid_.face_area(iface);

        // - 1/Re * Laplace(U)
        for (const std::pair<const size_t, double>& iter: df_dn_computer_.linear_combination(iface)) {
            size_t column = iter.first;
            double coef = -coef_diff * area * iter.second;
            mat.add_value(negative_colloc, column, coef);
            mat.add_value(positive_colloc, column, -coef);
        }

        // advection_u * grad u (using symmetric scheme)
        {
            double coef = area * un_face_[iface] / 2.0;
            if (!collocations_.is_boundary_colloc(negative_colloc)) {
                mat.add_value(negative_colloc, negative_colloc, coef);
                mat.add_value(negative_colloc, positive_colloc, coef);
            }
            if (!collocations_.is_boundary_colloc(positive_colloc)) {
                mat.add_value(positive_colloc, positive_colloc, -coef);
                mat.add_value(positive_colloc, negative_colloc, -coef);
            }
        }
    }
    return mat.to_csr();
}

void PoiseuilleFvmWorker::prepare_slae(const CsrMatrix& A) {
    mat_u_ = A;
    mat_v_ = A;
    double diag_relax = 1.0 / alpha_u_;
    double rhs_relax = (1.0 - alpha_u_) / alpha_u_;

    // internals
    for (auto [i, au_i, av_i]: matrix_iter::iv_diag(mat_u_, mat_v_)) {
        if (i < grid_.n_cells()) {
            // pressure equation coefficient
            du_[i] = grid_.cell_volume(i) * alpha_u_ / au_i;
            dv_[i] = grid_.cell_volume(i) * alpha_u_ / av_i;
            if (std::abs(du_[i] - dv_[i]) > 1e-6) {
                // dn computation is not ready for anisotropic diffusion coefficients
                _THROW_NOT_IMP_;
            }
            // rhs
            rhs_u_[i] = grid_.cell_volume(i) * (-grad_p_[i].x) + rhs_relax * au_i * u_[i];
            rhs_v_[i] = grid_.cell_volume(i) * (-grad_p_[i].y) + rhs_relax * av_i * v_[i];
            // matrix diag
            au_i *= diag_relax;
            av_i *= diag_relax;
        }
    }
    dn_face_ = compute_dn_face(du_, dv_);

    // Dirichlet boundary values
    // input bc: u=1, v=0
    for (size_t icolloc: boundary_info_.input) {
        mat_u_.set_unit_row(icolloc);
        rhs_u_[icolloc] = input_profile(collocations_.points[icolloc].y);
        mat_v_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
    }
    // wall bc: u=0
    for (size_t icolloc: boundary_info_.wall) {
        mat_u_.set_unit_row(icolloc);
        rhs_u_[icolloc] = 0.0;
        mat_v_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
    }
}

std::vector<double> PoiseuilleFvmWorker::compute_p(const std::vector<double>& un_face) const {
    LodMatrix Ap(vec_size());
    std::vector<double> rhs(vec_size(), 0);
    // assemble slae
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [negative_colloc, positive_colloc] = collocations_.tab_face_colloc(iface);
        double area = grid_.face_area(iface);
        // matrix: div(d*grad(p))
        for (const std::pair<const size_t, double>& iter: df_dn_computer_.linear_combination(iface)) {
            size_t column = iter.first;
            double coef = -dn_face_[iface] * area * iter.second;
            Ap.add_value(negative_colloc, column, coef);
            Ap.add_value(positive_colloc, column, -coef);
        }
        // rhs: -div(u)
        {
            double coef = area * un_face[iface];
            if (!collocations_.is_boundary_colloc(negative_colloc)) {
                rhs[negative_colloc] -= coef;
            }
            if (!collocations_.is_boundary_colloc(positive_colloc)) {
                rhs[positive_colloc] += coef;
            }
        }
    }
    // p = 0 at output
    for (size_t icolloc: boundary_info_.output) {
        Ap.set_unit_row(icolloc);
        rhs[icolloc] = 0;
    }
    // solve
    std::vector<double> p = p_;
    AmgcMatrixSolver::solve_slae(Ap.to_csr(), rhs, p);
    return p;
}

} // namespace

TEST_CASE("Poiseuille 2D, FVM-SIMPLE algorithm", "[pois-fvm]") {
    std::cout << std::endl << "--- cfd_test [pois-fvm] --- " << std::endl;

    // problem parameters
    double Re = 100;
    size_t max_it = 10000;
    double eps = 1e-1;
    double alpha_p = 0.3;
    double alpha_u = 0.8;

    // worker initialization
    size_t n_part = 20;
    RegularGrid2D grid(0, 2, -0.5, 0.5, 2 * n_part, n_part);
    PoiseuilleFvmWorker worker(grid, Re, alpha_u, alpha_p);
    worker.initialize_saver("poiseuille-fvm", 5);
    worker.save_current_fields(0);

    size_t it = 0;
    for (it = 1; it < max_it; ++it) {
        double nrm = worker.step();
        worker.save_current_fields(it);
        std::cout << it << ": " << nrm << std::endl;
        // break inner iterations if residual is low enough
        if (nrm < eps) {
            break;
        } else if (it == max_it - 1) {
            std::cout << "WARNING: internal SIMPLE interations not converged with nrm = " << nrm << std::endl;
        }
    }
    CHECK(it == 19);
}
