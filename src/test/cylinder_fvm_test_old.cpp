#include "cfd/debug/printer.hpp"
#include "cfd/debug/saver.hpp"
#include "cfd/debug/tictoc.hpp"
#include "cfd/fvm/fvm_dfdn.hpp"
#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/fvm/fvm_gradient.hpp"
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
    CylinderFvmWorker(const IGrid& grid, double Re, double alpha_u, double alpha_p);
    void initialize_saver(std::string stem, double save_time_step);

    double step();
    void save_current_fields(double time) const;
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
        std::vector<size_t> cyl;
        std::vector<size_t> input;
        std::vector<size_t> output;
        std::vector<size_t> sym;
    };
    BoundaryInfo boundary_info_;

    std::vector<double> p_;
    std::vector<double> u_;
    std::vector<double> v_;

    CsrMatrix mat_u_;
    CsrMatrix mat_v_;
    std::vector<double> rhs_u_;
    std::vector<double> rhs_v_;
    std::vector<double> du_;
    std::vector<double> dv_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_;

    CsrMatrix assemble_u_matrix() const;
    void prepare_slae(const CsrMatrix& A);
    std::vector<double> compute_u_tilde(const CsrMatrix& A, const std::vector<double>& d,
                                        const std::vector<double>& u_star, const std::vector<double>& u) const;
    std::vector<double> compute_p_tilde(const std::vector<double>& u_tilde, const std::vector<double>& v_tilde) const;
    std::vector<double> compute_un_face(const std::vector<double>& x, const std::vector<double>& y) const;
    std::vector<double> compute_dn_face(const std::vector<double>& x, const std::vector<double>& y) const;

    double to_next_iteration();
    void gather_boundary_collocations();
};

CylinderFvmWorker::CylinderFvmWorker(const IGrid& grid, double Re, double alpha_u, double alpha_p)
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

    to_next_iteration();
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

double CylinderFvmWorker::to_next_iteration() {
    CsrMatrix A = assemble_u_matrix();
    prepare_slae(A);

    // residual vectors
    std::vector<double> res_u = compute_residual_vec(mat_u_, rhs_u_, u_);
    std::vector<double> res_v = compute_residual_vec(mat_v_, rhs_v_, v_);
    // norm
    double res = 0;
    for (size_t icell = 0; icell < grid_.n_cells(); ++icell) {
        res = std::max(res, std::max(res_u[icell], res_v[icell]));
    }
    return res;
};

double CylinderFvmWorker::step() {
    // 1. Predictor step: U-star
    std::vector<double> u_star = u_, v_star = v_;
    AmgcMatrixSolver::solve_slae(mat_u_, rhs_u_, u_star);
    AmgcMatrixSolver::solve_slae(mat_v_, rhs_v_, v_star);
    // 2. Rhi-Chow velocity
    std::vector<double> u_tilde = compute_u_tilde(mat_u_, du_, u_star, u_);
    std::vector<double> v_tilde = compute_u_tilde(mat_v_, dv_, v_star, v_);
    // 3. New pressure
    std::vector<double> p_tilde = compute_p_tilde(u_tilde, v_tilde);
    // 4. New velocity
    // -- internal
    for (size_t icell = 0; icell < grid_.n_cells(); ++icell) {
        Vector grad_p_tilde = grad_computer_->cell_compute(icell, p_tilde);
        u_[icell] = u_tilde[icell] - du_[icell] * grad_p_tilde.x;
        v_[icell] = v_tilde[icell] - dv_[icell] * grad_p_tilde.y;
    }
    // -- boundary values
    for (size_t i = grid_.n_cells(); i < vec_size(); ++i) {
        u_[i] = u_tilde[i];
        v_[i] = v_tilde[i];
    }
    // 5. Relax Pressure
    p_ = vector_sum(1 - alpha_p_, p_, alpha_p_, p_tilde);

    return to_next_iteration();
}

void CylinderFvmWorker::save_current_fields(double time) const {
    if (writer_) {
        std::string filepath = writer_->add(time);
        if (filepath.empty()) {
            return;
        }
        grid_.save_vtk(filepath);
        VtkUtils::add_cell_data(p_, "pressure", filepath, grid_.n_cells());
        VtkUtils::add_cell_vector(u_, v_, "velocity", filepath, grid_.n_cells());
    }
}

std::vector<double> CylinderFvmWorker::compute_un_face(const std::vector<double>& x,
                                                       const std::vector<double>& y) const {
    std::vector<double> ret(grid_.n_faces(), 0.0);
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [ci, cj] = collocations_.tab_face_colloc(iface);
        Vector normal = grid_.face_normal(iface);
        Vector u_face;
        if (collocations_.is_boundary_colloc(ci)) {
            u_face = Vector{x[ci], y[ci]};
        } else if (collocations_.is_boundary_colloc(cj)) {
            u_face = Vector{x[cj], y[cj]};
        } else {
            u_face = Vector{0.5 * (x[ci] + x[cj]), 0.5 * (y[ci] + y[cj])};
        }
        ret[iface] = dot_product(u_face, normal);
    }
    return ret;
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
        ret[iface] = dot_product(u_face, normal);
    }
    return ret;
}

CsrMatrix CylinderFvmWorker::assemble_u_matrix() const {
    LodMatrix mat(vec_size());
    std::vector<double> un_face = compute_un_face(u_, v_);
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
            double coef = area * un_face[iface] / 2.0;
            mat.add_value(negative_colloc, negative_colloc, coef);
            mat.add_value(negative_colloc, positive_colloc, coef);
            mat.add_value(positive_colloc, positive_colloc, -coef);
            mat.add_value(positive_colloc, negative_colloc, -coef);
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
            du_[i] = alpha_u_ / au_i;
            dv_[i] = alpha_u_ / av_i;
            // rhs
            Vector grad_p = grad_computer_->cell_compute(i, p_);
            rhs_u_[i] = grid_.cell_volume(i) * (-grad_p.x + rhs_relax * au_i * u_[i]);
            rhs_v_[i] = grid_.cell_volume(i) * (-grad_p.y + rhs_relax * av_i * v_[i]);
            // matrix diag
            au_i *= diag_relax;
            av_i *= diag_relax;
        }
    }

    // Dirichlet boundary values
    // input bc: u=1, v=0
    for (size_t icolloc: boundary_info_.input) {
        mat_u_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
        mat_v_.set_unit_row(icolloc);
        rhs_u_[icolloc] = 1.0;
    }
    // sym: v = 0
    for (size_t icolloc: boundary_info_.sym) {
        mat_u_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
    }
    // cyl: u = v = 0
    for (size_t icolloc: boundary_info_.cyl) {
        mat_u_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
        mat_v_.set_unit_row(icolloc);
        rhs_u_[icolloc] = 0.0;
    }
    // output bc: u=1, v=0;
    for (size_t icolloc: boundary_info_.output) {
        mat_u_.set_unit_row(icolloc);
        rhs_v_[icolloc] = 0.0;
        mat_v_.set_unit_row(icolloc);
        rhs_u_[icolloc] = 0.0;
    }
}

std::vector<double> CylinderFvmWorker::compute_u_tilde(const CsrMatrix& A, const std::vector<double>& d,
                                                       const std::vector<double>& u_star,
                                                       const std::vector<double>& u) const {
    std::vector<double> u_tilde(vec_size());
    // internal
    for (size_t i = 0; i < grid_.n_cells(); ++i) {
        u_tilde[i] = d[i] * (-nodiag_mult(i, A, u_star)) + (1 - alpha_u_) * u[i];
    }
    // boundary
    for (size_t icolloc: boundary_info_.all) {
        u_tilde[icolloc] = u_star[icolloc];
    }
    return u_tilde;
}

std::vector<double> CylinderFvmWorker::compute_p_tilde(const std::vector<double>& u_tilde,
                                                       const std::vector<double>& v_tilde) const {
    LodMatrix Ap(vec_size());
    std::vector<double> rhs(vec_size(), 0);
    std::vector<double> dn_face = compute_dn_face(du_, dv_);
    std::vector<double> un_tilde_face = compute_un_face(u_tilde, v_tilde);
    // assemble slae
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [negative_colloc, positive_colloc] = collocations_.tab_face_colloc(iface);
        double area = grid_.face_area(iface);
        // matrix: div(d*grad(p))
        for (const std::pair<const size_t, double>& iter: df_dn_computer_.linear_combination(iface)) {
            size_t column = iter.first;
            double coef = -dn_face[iface] * area * iter.second;
            Ap.add_value(negative_colloc, column, coef);
            Ap.add_value(positive_colloc, column, -coef);
        }
        // rhs: -div(u_tilde)
        {
            double coef = area * un_tilde_face[iface];
            if (!collocations_.is_boundary_colloc(negative_colloc)) {
                rhs[negative_colloc] -= coef;
            }
            if (!collocations_.is_boundary_colloc(positive_colloc)) {
                rhs[positive_colloc] += coef;
            }
        }
    }
    Ap.set_unit_row(0);
    rhs[0] = 0.0;
    // solve
    std::vector<double> p_tilde = p_;
    AmgcMatrixSolver::solve_slae(Ap.to_csr(), rhs, p_tilde);
    return p_tilde;
}

std::string convergence_report(double time, size_t it) {
    std::ostringstream oss;
    oss << std::setprecision(2) << std::fixed;
    oss << "Time = " << std::setw(5) << time << " converged in " << it << " iterations" << std::endl;
    return oss.str();
}

} // namespace

TEST_CASE("Cylinder, FVM-SIMPLE algorithm", "[cylinder-fvm]") {
    std::cout << std::endl << "--- cfd_test [cylinder-fvm] --- " << std::endl;

    // problem parameters
    double Re = 30;
    size_t max_it = 10000;
    double alpha_p = 0.3;
    double alpha_u = 1.0;
    double time_step = 0.5;
    double end_time = 0.5;
    double eps = 0.1;

    // worker initialization
    auto grid = UnstructuredGrid2D::vtk_read(test_directory_file("cylgrid_5k.vtk"));
    CylinderFvmWorker worker(grid, Re, alpha_u, alpha_p);
    worker.initialize_saver("cylinder2-fvm", 0.5);
    worker.save_current_fields(0.0);

    // time loop
    size_t it = 0;
    for (double time = time_step; time < end_time + 1e-6; time += time_step) {
        for (it = 1; it < max_it; ++it) {
            double nrm = worker.step();
            // break inner iterations if residual is low enough
            if (nrm < eps) {
                break;
            } else if (it == max_it - 1) {
                std::cout << "WARNING: internal SIMPLE interations not converged with nrm = " << nrm << std::endl;
            }
            std::cout << it << ": " << nrm << std::endl;
            worker.save_current_fields(static_cast<double>(it));
        }
        // uvp_old = uvp
        worker.to_next_time_step();

        // save and report
        std::cout << convergence_report(time, it);
        worker.save_current_fields(time);
    }

    CHECK(it == 13);
}
