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

// => (M - diag(M))* u
double nodiag_mult(size_t irow, const CsrMatrix& mat, const std::vector<double>& u) {
    double ret = 0;
    for (size_t a = mat.addr()[irow]; a < mat.addr()[irow + 1]; ++a) {
        if (mat.cols()[a] != irow) {
            ret += mat.vals()[a] * u[mat.cols()[a]];
        }
    }
    return ret;
}

struct CylinderFvmWorker {
    CylinderFvmWorker(const IGrid& grid, double Re, double time_step, double alpha_u = 1.0, double alpha_p = 1.0);
    void initialize_saver(std::string stem, double save_time_step);

    void step(size_t n_piso, double piso_epsilon);
    void save_current_fields(double time, bool force = false) const;
    void init_time_step();
    double init_step();

    size_t vec_size() const {
        return collocations_.size();
    }

    double testing_value() const {
        return vector_abs(Vector{u_[954], v_[954]});
    }

private:
    const IGrid& grid_;
    const double Re_;
    const double dt_;
    const double alpha_u_;
    const double alpha_p_;
    const FvmExtendedCollocations collocations_;
    const FvmFacesDn df_dn_computer_;
    const std::shared_ptr<IFvmCellGradient> grad_computer_;

    struct BoundaryInfo {
        std::vector<size_t> all;
        std::vector<size_t> sym;
        std::vector<size_t> input;
        std::vector<size_t> output;
        std::vector<size_t> cyl;
    };
    BoundaryInfo boundary_info_;

    std::vector<double> p_;
    std::vector<double> u_, v_;
    std::vector<double> u_old_, v_old_;
    std::vector<double> du_, dv_;
    std::vector<double> un_face_;
    std::vector<double> dn_face_;
    std::vector<double> dpdn_face_;
    std::vector<Vector> grad_p_;

    CsrMatrix mat_u_, mat_v_;
    std::vector<double> rhs_u_, rhs_v_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_;

    CsrMatrix assemble_u_matrix() const;
    AmgcMatrixSolver assemble_p_solver() const;
    void prepare_slae(const CsrMatrix& A);
    std::vector<double> compute_u_tilde(const CsrMatrix& A, const std::vector<double>& d, double time_step,
                                        const std::vector<double>& u_star, const std::vector<double>& u) const;
    std::vector<double> compute_p(const AmgcMatrixSolver& solver, const std::vector<double>& un_face) const;
    std::vector<double> compute_dn_face(const std::vector<double>& x, const std::vector<double>& y) const;
    std::vector<double> compute_un_face(const std::vector<double>& u, const std::vector<double>& v) const;
    std::vector<double> compute_un_face_rhie_chow(const std::vector<double>& u_star,
                                                  const std::vector<double>& v_star) const;

    void gather_boundary_collocations();
};

CylinderFvmWorker::CylinderFvmWorker(const IGrid& grid, double Re, double time_step, double alpha_u, double alpha_p)
    : grid_(grid),
      Re_(Re),
      dt_(time_step),
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
    u_old_ = std::vector<double>(vec_size(), 1);
    v_old_ = std::vector<double>(vec_size(), 0);
    p_ = std::vector<double>(vec_size(), 0);
    un_face_ = std::vector<double>(grid_.n_faces(), 0);

    dn_face_ = std::vector<double>(grid_.n_faces(), 0);
    dpdn_face_ = std::vector<double>(grid_.n_faces(), 0);
    grad_p_ = std::vector<Vector>(grid_.n_cells(), Vector{0, 0, 0});
}

void CylinderFvmWorker::gather_boundary_collocations() {
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
            bi.sym.push_back(icolloc);
        } else if (std::abs(fc.y - ymax) < 1e-6) {
            bi.sym.push_back(icolloc);
        } else if (std::abs(fc.x - xmin) < 1e-6) {
            bi.input.push_back(icolloc);
        } else if (std::abs(fc.x - xmax) < 1e-6) {
            bi.output.push_back(icolloc);
        } else {
            bi.cyl.push_back(icolloc);
        }
    }
}

void CylinderFvmWorker::initialize_saver(std::string stem, double save_time_step) {
    writer_.reset(new VtkUtils::TimeSeriesWriter(stem));
    writer_->set_time_step(save_time_step);
};

void CylinderFvmWorker::init_time_step() {
    u_old_ = u_;
    v_old_ = v_;
}

double CylinderFvmWorker::init_step() {
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

void CylinderFvmWorker::step(size_t n_piso, double piso_epsilon) {
    // 1. Predictor step: U-star
    std::vector<double> u_star = u_, v_star = v_;
    AmgcMatrixSolver::solve_slae(mat_u_, rhs_u_, u_star);
    AmgcMatrixSolver::solve_slae(mat_v_, rhs_v_, v_star);
    std::vector<double> un_star_face = compute_un_face_rhie_chow(u_star, v_star);

    // 2. Pressure correction
    AmgcMatrixSolver p_solver = assemble_p_solver();
    std::vector<double> p_prime = compute_p(p_solver, un_star_face);
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

    // 4. Set values for SIMPLE step
    std::vector<double> u_new = vector_sum(u_star, u_prime);
    std::vector<double> v_new = vector_sum(v_star, v_prime);
    std::vector<double> un_face_new = vector_sum(un_star_face, un_prime_face);
    std::vector<double> p_new = vector_sum(p_, alpha_p_, p_prime);
    std::vector<Vector> grad_p_new = vector_sum(grad_p_, alpha_p_, grad_p_prime);
    std::vector<double> dpdn_face_new = vector_sum(dpdn_face_, alpha_p_, dpprime_dn_face);

    // 5. PISO steps
    for (size_t ipiso = 0; true; ++ipiso) {
        // compute piso error
        double res = 0.0;
        for (size_t icell = 0; icell < grid_.n_cells(); ++icell) {
            Vector grad = grad_p_new[icell] - grad_computer_->cell_compute(icell, p_);
            double err_u =
                mat_u_.mult_vec(icell, u_new) - mat_u_.mult_vec(icell, u_star) + grid_.cell_volume(icell) * grad.x;
            double err_v =
                mat_v_.mult_vec(icell, v_new) - mat_v_.mult_vec(icell, v_star) + grid_.cell_volume(icell) * grad.y;
            res = std::max(res, std::max(std::abs(err_u), std::abs(err_v)));
        }
        std::cout << "--- " << ipiso << ": PISO error = " << res << " iterations." << std::endl;
        if (res < piso_epsilon) {
            std::cout << "--- PISO converged in " << ipiso << " iterations." << std::endl;
            break;
        } else if (ipiso >= n_piso) {
            std::cout << "--- WARNING: PISO failed to converge with residual = " << res << std::endl;
            break;
        }
        // u_tilde
        std::vector<double> u_tilde(vec_size(), 0), v_tilde(vec_size(), 0);
        for (size_t i = 0; i < grid_.n_cells(); ++i) {
            u_tilde[i] = -du_[i] / grid_.cell_volume(i) * nodiag_mult(i, mat_u_, u_prime);
            v_tilde[i] = -dv_[i] / grid_.cell_volume(i) * nodiag_mult(i, mat_v_, v_prime);
        }
        std::vector<double> un_tilde_face = compute_un_face(u_tilde, v_tilde);
        // p_prime
        p_prime = compute_p(p_solver, un_tilde_face);
        grad_p_prime = grad_computer_->compute(p_prime);
        dpprime_dn_face = df_dn_computer_.compute(p_prime);
        // u_prime
        for (size_t i = 0; i < grid_.n_cells(); ++i) {
            u_prime[i] = u_tilde[i] - du_[i] * grad_p_prime[i].x;
            v_prime[i] = v_tilde[i] - dv_[i] * grad_p_prime[i].y;
        }
        for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
            un_prime_face[iface] = un_tilde_face[iface] - dn_face_[iface] * dpprime_dn_face[iface];
        }
        // add to new values
        u_new = vector_sum(u_new, u_prime);
        v_new = vector_sum(v_new, v_prime);
        un_face_new = vector_sum(un_face_new, un_prime_face);
        p_new = vector_sum(p_new, alpha_p_, p_prime);
        grad_p_new = vector_sum(grad_p_new, alpha_p_, grad_p_prime);
        dpdn_face_new = vector_sum(dpdn_face_new, alpha_p_, dpprime_dn_face);
    }

    // 6. Set final values
    std::swap(u_, u_new);
    std::swap(v_, v_new);
    std::swap(un_face_, un_face_new);
    std::swap(p_, p_new);
    std::swap(grad_p_, grad_p_new);
    std::swap(dpdn_face_, dpdn_face_new);
}

void CylinderFvmWorker::save_current_fields(double time, bool force) const {
    if (writer_) {
        std::string filepath = writer_->add(time, force);
        if (filepath.empty()) {
            return;
        }
        grid_.save_vtk(filepath);
        VtkUtils::add_cell_data(p_, "pressure", filepath, grid_.n_cells());
        VtkUtils::add_cell_vector(u_, v_, "velocity", filepath, grid_.n_cells());
    }
}

std::vector<double> CylinderFvmWorker::compute_dn_face(const std::vector<double>& x,
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
        // d is a tensor tensor(x,0; 0, y) => d_nn = d_xx * nx^2 + d_yy * ny^2, d_ns = 0 assuming d_xx = d_yy
        ret[iface] = u_face.x * normal.x * normal.x + u_face.y * normal.y * normal.y;
    }
    return ret;
}

std::vector<double> CylinderFvmWorker::compute_un_face(const std::vector<double>& u,
                                                       const std::vector<double>& v) const {

    std::vector<double> ret(grid_.n_faces());

    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [ci, cj] = grid_.tab_face_cell(iface);
        if (ci == INVALID_INDEX || cj == INVALID_INDEX) {
            size_t icolloc = collocations_.boundary_colloc(iface);
            ret[iface] = dot_product(Vector{u[icolloc], v[icolloc]}, grid_.face_normal(iface));
        } else {
            Vector normal = grid_.face_normal(iface);
            Vector uvec_i = Vector(u[ci], v[ci]);
            Vector uvec_j = Vector(u[cj], v[cj]);

            double u_i = dot_product(uvec_i, normal);
            double u_j = dot_product(uvec_j, normal);

            ret[iface] = 0.5 * (u_i + u_j);
        }
    }

    return ret;
}

std::vector<double> CylinderFvmWorker::compute_un_face_rhie_chow(const std::vector<double>& u,
                                                                 const std::vector<double>& v) const {

    std::vector<double> ret(grid_.n_faces());

    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [ci, cj] = grid_.tab_face_cell(iface);
        if (ci == INVALID_INDEX || cj == INVALID_INDEX) {
            size_t icolloc = collocations_.boundary_colloc(iface);
            ret[iface] = dot_product(Vector{u[icolloc], v[icolloc]}, grid_.face_normal(iface));
        } else {
            // Rhie-Chow interpolation
            Vector normal = grid_.face_normal(iface);
            Vector uvec_i = Vector(u[ci], v[ci]);
            Vector uvec_j = Vector(u[cj], v[cj]);

            double u_i = dot_product(uvec_i, normal);
            double u_j = dot_product(uvec_j, normal);
            double d_dpdn_i = dot_product(Vector{grad_p_[ci].x * du_[ci], grad_p_[ci].y * dv_[ci]}, normal);
            double d_dpdn_j = dot_product(Vector{grad_p_[cj].x * du_[cj], grad_p_[cj].y * dv_[cj]}, normal);
            double d_dpdn_ij = dn_face_[iface] * dpdn_face_[iface];

            ret[iface] = 0.5 * (u_i + u_j) + 0.5 * (d_dpdn_i + d_dpdn_j - 2 * d_dpdn_ij);
        }
    }

    return ret;
}

CsrMatrix CylinderFvmWorker::assemble_u_matrix() const {
    LodMatrix mat(vec_size());
    double coef_diff = 1.0 / Re_;

    // du / dt only for internal collocations
    for (size_t icell = 0; icell < grid_.n_cells(); ++icell) {
        mat.add_value(icell, icell, grid_.cell_volume(icell) / dt_);
    }
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [negative_colloc, positive_colloc] = collocations_.tab_face_colloc(iface);
        double area = grid_.face_area(iface);

        // - 1/Re * Laplace(U) including (-1/Re du/dn) approx on boundary collocations
        for (const std::pair<const size_t, double>& iter: df_dn_computer_.linear_combination(iface)) {
            size_t column = iter.first;
            double coef = -coef_diff * area * iter.second;
            mat.add_value(negative_colloc, column, coef);
            mat.add_value(positive_colloc, column, -coef);
        }

        // advection_u * grad u (using symmetric scheme) only for internal collocations
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

void CylinderFvmWorker::prepare_slae(const CsrMatrix& A) {
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
            rhs_u_[i] = grid_.cell_volume(i) * (u_old_[i] / dt_ - grad_p_[i].x) + rhs_relax * au_i * u_[i];
            rhs_v_[i] = grid_.cell_volume(i) * (v_old_[i] / dt_ - grad_p_[i].y) + rhs_relax * av_i * v_[i];
            // matrix diag relax
            au_i *= diag_relax;
            av_i *= diag_relax;
        }
    }
    dn_face_ = compute_dn_face(du_, dv_);

    // Dirichlet boundary values
    // input bc: u=1, v=0
    for (size_t icolloc: boundary_info_.input) {
        mat_u_.set_unit_row(icolloc);
        rhs_u_[icolloc] = 1.0;
        mat_v_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
    }
    // sym bc: v=0
    for (size_t icolloc: boundary_info_.sym) {
        mat_v_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
    }
    // cyl bc: u=v=0
    for (size_t icolloc: boundary_info_.cyl) {
        mat_u_.set_unit_row(icolloc);
        rhs_u_[icolloc] = 0.0;
        mat_v_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
    }
}

AmgcMatrixSolver CylinderFvmWorker::assemble_p_solver() const {
    LodMatrix Ap(vec_size());
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
    }
    // p = 0 at output
    for (size_t icolloc: boundary_info_.output) {
        Ap.set_unit_row(icolloc);
    }
    return AmgcMatrixSolver(Ap.to_csr());
}

std::vector<double> CylinderFvmWorker::compute_p(const AmgcMatrixSolver& solver,
                                                 const std::vector<double>& un_face) const {
    std::vector<double> rhs(vec_size(), 0);
    // assemble rhs = -div(u)
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [negative_colloc, positive_colloc] = collocations_.tab_face_colloc(iface);
        double coef = grid_.face_area(iface) * un_face[iface];
        if (!collocations_.is_boundary_colloc(negative_colloc)) {
            rhs[negative_colloc] -= coef;
        }
        if (!collocations_.is_boundary_colloc(positive_colloc)) {
            rhs[positive_colloc] += coef;
        }
    }
    // p = 0 at output
    for (size_t icolloc: boundary_info_.output) {
        rhs[icolloc] = 0;
    }
    // solve
    std::vector<double> p;
    solver.solve(rhs, p);
    return p;
}

} // namespace

TEST_CASE("Cylinder 2D, FVM-SIMPLE algorithm", "[cylinder-fvm]") {
    std::cout << std::endl << "--- cfd_test [cylinder-fvm] --- " << std::endl;

    // problem parameters
    double Re = 100;
    double time_step = 0.05;
    double end_time = 2.0;
    size_t max_piso_iters = 10;
    double piso_eps = 1e-3;

    // worker initialization
    auto grid = UnstructuredGrid2D::vtk_read(test_directory_file("cylgrid_5k.vtk"));
    CylinderFvmWorker worker(grid, Re, time_step);
    worker.initialize_saver("cylinder-fvm", 0.5);
    worker.save_current_fields(0);

    double time = 0.0;
    for (time = time_step; time < end_time + 1e-6; time += time_step) {
        std::cout << "time= " << time << std::endl;
        worker.init_time_step();
        worker.init_step();
        worker.step(max_piso_iters, piso_eps);
        worker.save_current_fields(time);
    }
    worker.save_current_fields(time, true);

    CHECK(worker.testing_value() == Approx(1.1191340904).epsilon(1e-3));
}
