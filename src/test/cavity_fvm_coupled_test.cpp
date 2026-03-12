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
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"
#include "utils/filesystem.hpp"
#include "utils/vecmat.hpp"
#include <list>

using namespace cfd;

namespace {

double MC(double r) {
    return std::max(0.0, std::min(std::min(2.0 * r, (1 + r) / 2.0), 2.0));
}

double tvd_upwind_weight(const double* u, const std::vector<Vector>& grad_u, Vector cij, size_t i, size_t j) {

    // treat boundary: take boundary value
    size_t n_cells = grad_u.size();
    if (i >= n_cells) {
        return 1.0;
    } else if (j >= n_cells) {
        return 0.0;
    }

    // internals
    double dudc = dot_product(grad_u[i], cij);
    double up = u[j] - 2 * dudc;
    double denum = u[i] - u[j];
    if (std::abs(denum) < 1e-12)
        denum = 1e-12;
    double r = (u[i] - up) / denum;
    double phi = MC(r);
    return 1 - phi / 2;
}

} // namespace

struct CavityFvmCoupledWorker {
    CavityFvmCoupledWorker(const IGrid& grid, double Re, double tau);
    void initialize_saver(std::string stem);
    double step();
    void save_current_fields(size_t iter);

    size_t vec_size() const {
        return collocations_.size();
    }

private:
    const IGrid& grid_;
    const double Re_;
    const double delta_t_;
    const FvmExtendedCollocations collocations_;
    const FvmFacesDn dfdn_computer_;
    const GaussLinearFvmCellGradient grad_computer_;

    struct BoundaryInfo {
        std::vector<size_t> bnd_colloc;
        std::vector<size_t> bnd_colloc_u0;
        std::vector<size_t> bnd_colloc_u1;
    };
    BoundaryInfo boundary_info_;

    std::vector<double> un_face_;
    std::vector<Vector> grad_p_;
    std::vector<double> du_, dv_, dn_face_;

    CsrMatrix mat_;
    std::vector<double> rhs_;
    std::vector<double> uvp_;

    double* u_;
    double* v_;
    double* p_;

    std::shared_ptr<VtkUtils::TimeSeriesWriter> writer_;

    void to_next_iteration();
    void gather_boundary_collocations();
    std::vector<double> compute_un_rhie_chow() const;
    void compute_d(const CsrMatrix& Bu, const CsrMatrix& Bv);
    void u_equation(CsrMatrix&, std::vector<double>&) const;
    void v_equation(CsrMatrix&, std::vector<double>&) const;
    void continuity_equation(CsrMatrix&, std::vector<double>&) const;
    void assemble_slae();
};

CavityFvmCoupledWorker::CavityFvmCoupledWorker(const IGrid& grid, double Re, double tau)
    : grid_(grid),
      Re_(Re),
      delta_t_(tau),
      collocations_(grid),
      dfdn_computer_(grid, collocations_),
      grad_computer_(grid, collocations_) {
    gather_boundary_collocations();

    uvp_ = std::vector<double>(3 * vec_size(), 0);
    rhs_ = std::vector<double>(3 * vec_size(), 0);
    u_ = uvp_.data();
    v_ = uvp_.data() + vec_size();
    p_ = uvp_.data() + 2 * vec_size();
    un_face_ = std::vector<double>(grid_.n_faces(), 0);
    grad_p_ = std::vector<Vector>(grid_.n_cells(), {0, 0, 0});
    du_ = std::vector<double>(vec_size(), 0);
    dv_ = std::vector<double>(vec_size(), 0);
    dn_face_ = std::vector<double>(grid_.n_faces(), 0);
    to_next_iteration();
}

void CavityFvmCoupledWorker::gather_boundary_collocations() {
    BoundaryInfo& bi = boundary_info_;

    std::list<std::pair<size_t, size_t>> colloc_faces;
    for (size_t icolloc: collocations_.face_collocations) {
        size_t iface = collocations_.face_index(icolloc);
        bi.bnd_colloc.push_back(icolloc);
        if (std::abs(grid_.face_center(iface).y - 1) < 1e-6) {
            bi.bnd_colloc_u1.push_back(icolloc);
        } else {
            bi.bnd_colloc_u0.push_back(icolloc);
        }
    }
}

void CavityFvmCoupledWorker::initialize_saver(std::string stem) {
    writer_.reset(new VtkUtils::TimeSeriesWriter(stem));
};

void CavityFvmCoupledWorker::to_next_iteration() {
    assemble_slae();
};

void CavityFvmCoupledWorker::compute_d([[maybe_unused]] const CsrMatrix& Bu, [[maybe_unused]] const CsrMatrix& Bv) {
    // du, dv
    for (size_t i = 0; i < grid_.n_cells(); ++i) {
        double vol = grid_.cell_volume(i);
        du_[i] = vol / Bu.value(i, i);
        dv_[i] = vol / Bv.value(i, vec_size() + i);
    }

    // dn_face
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        auto [ci, cj] = collocations_.tab_face_colloc(iface);
        Vector normal = grid_.face_normal(iface);
        Vector d_face;
        if (collocations_.is_boundary_colloc(ci)) {
            d_face = Vector{du_[cj], dv_[cj]};
        } else if (collocations_.is_boundary_colloc(cj)) {
            d_face = Vector{du_[ci], dv_[ci]};
        } else {
            d_face = Vector{0.5 * (du_[ci] + du_[cj]), 0.5 * (dv_[ci] + dv_[cj])};
        }
        // d is a tensor tensor(x,0; 0, y) => d_nn = d_xx * nx^2 + d_yy * ny^2, d_ns = 0 assuming d_xx = d_yy
        dn_face_[iface] = d_face.x * normal.x * normal.x + d_face.y * normal.y * normal.y;
    }
    du_.back() = 0;
}

void CavityFvmCoupledWorker::u_equation(CsrMatrix& B, std::vector<double>& rhs) const {
    LodMatrix B0(vec_size()), B2(vec_size());
    std::fill(rhs.begin(), rhs.end(), 0.0);
    std::vector<Vector> grad_u = grad_computer_.compute(u_);

    // cell loop
    for (size_t icell = 0; icell < grid_.n_cells(); ++icell) {
        double v = grid_.cell_volume(icell) / delta_t_;
        // nonstat
        B0.add_value(icell, icell, v);
        rhs[icell] += v * u_[icell];
    }
    // face loop
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        size_t negative_colloc = collocations_.tab_face_colloc(iface)[0];
        size_t positive_colloc = collocations_.tab_face_colloc(iface)[1];
        double area = grid_.face_area(iface);
        double un = un_face_[iface];

        // - 1/Re * Laplace(u)
        for (const std::pair<const size_t, double>& iter: dfdn_computer_.linear_combination(iface)) {
            size_t column = iter.first;
            double coef = -1.0 / Re_ * area * iter.second;
            B0.add_value(negative_colloc, column, coef);
            B0.add_value(positive_colloc, column, -coef);
        }

        // + nonlinear tvd convection
        {
            size_t i = negative_colloc; // upwind cell
            size_t j = positive_colloc; // downwind cell
            double coef = area * un;
            if (un < 0) {
                std::swap(i, j);
                coef *= -1;
            };
            Vector cij = collocations_.points[j] - collocations_.points[i];
            double wi = tvd_upwind_weight(u_, grad_u, cij, i, j);
            double wj = 1 - wi;
            B0.add_value(i, i, wi * coef);
            B0.add_value(i, j, wj * coef);
            B0.add_value(j, j, -wj * coef);
            B0.add_value(j, i, -wi * coef);
        }
        // + dp/dx
        {
            // area * nx / 2
            double value = area * grid_.face_normal(iface).x * 0.5;

            if (collocations_.is_boundary_colloc(negative_colloc)) {
                B2.add_value(positive_colloc, negative_colloc, -2 * value);
            } else if (collocations_.is_boundary_colloc(positive_colloc)) {
                B2.add_value(negative_colloc, positive_colloc, 2 * value);
            } else {
                B2.add_value(negative_colloc, negative_colloc, value);
                B2.add_value(negative_colloc, positive_colloc, value);
                B2.add_value(positive_colloc, negative_colloc, -value);
                B2.add_value(positive_colloc, positive_colloc, -value);
            }
        }
    }
    // boundary conditions
    for (size_t icolloc: boundary_info_.bnd_colloc) {
        B0.set_unit_row(icolloc);
        B2.remove_row(icolloc);
        rhs[icolloc] = 0;
    }
    for (size_t icolloc: boundary_info_.bnd_colloc_u1) {
        rhs[icolloc] = 1;
    }

    // block assemble
    B = assemble_block_matrix(vec_size(), vec_size(), {{&B0, nullptr, &B2}});
};

void CavityFvmCoupledWorker::v_equation(CsrMatrix& B, std::vector<double>& rhs) const {
    LodMatrix B1(vec_size()), B2(vec_size());
    std::fill(rhs.begin(), rhs.end(), 0.0);
    std::vector<Vector> grad_v = grad_computer_.compute(v_);

    // cell loop
    for (size_t icell = 0; icell < grid_.n_cells(); ++icell) {
        double v = grid_.cell_volume(icell) / delta_t_;
        // nonstat
        B1.add_value(icell, icell, v);
        rhs[icell] += v * v_[icell];
    }
    // face loop
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        size_t negative_colloc = collocations_.tab_face_colloc(iface)[0];
        size_t positive_colloc = collocations_.tab_face_colloc(iface)[1];
        double area = grid_.face_area(iface);
        double un = un_face_[iface];

        // - 1.0/Re * Laplace(u)
        for (const std::pair<const size_t, double>& iter: dfdn_computer_.linear_combination(iface)) {
            size_t column = iter.first;
            double coef = -1.0 / Re_ * area * iter.second;
            B1.add_value(negative_colloc, column, coef);
            B1.add_value(positive_colloc, column, -coef);
        }

        // + nonlinear tvd convection
        {
            size_t i = negative_colloc; // upwind cell
            size_t j = positive_colloc; // downwind cell
            double coef = area * un;
            if (un < 0) {
                std::swap(i, j);
                coef *= -1;
            };
            Vector cij = collocations_.points[j] - collocations_.points[i];
            double wi = tvd_upwind_weight(v_, grad_v, cij, i, j);
            double wj = 1 - wi;
            B1.add_value(i, i, wi * coef);
            B1.add_value(i, j, wj * coef);
            B1.add_value(j, j, -wj * coef);
            B1.add_value(j, i, -wi * coef);
        }
        // + dp/dy
        {
            // delta_t_ * area * ny / 2
            double value = area * grid_.face_normal(iface).y * 0.5;

            if (collocations_.is_boundary_colloc(negative_colloc)) {
                B2.add_value(positive_colloc, negative_colloc, -2 * value);
            } else if (collocations_.is_boundary_colloc(positive_colloc)) {
                B2.add_value(negative_colloc, positive_colloc, 2 * value);
            } else {
                B2.add_value(negative_colloc, negative_colloc, value);
                B2.add_value(negative_colloc, positive_colloc, value);
                B2.add_value(positive_colloc, negative_colloc, -value);
                B2.add_value(positive_colloc, positive_colloc, -value);
            }
        }
    }

    // boundary conditions
    for (size_t icolloc: boundary_info_.bnd_colloc) {
        B1.set_unit_row(icolloc);
        B2.remove_row(icolloc);
        rhs[icolloc] = 0;
    }

    // block assemble
    B = assemble_block_matrix(vec_size(), vec_size(), {{nullptr, &B1, &B2}});
};

void CavityFvmCoupledWorker::continuity_equation(CsrMatrix& B, std::vector<double>& rhs) const {
    LodMatrix B0(vec_size()), B1(vec_size()), B2(vec_size());
    std::fill(rhs.begin(), rhs.end(), 0.0);

    // face loop
    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        size_t negative_colloc = collocations_.tab_face_colloc(iface)[0];
        size_t positive_colloc = collocations_.tab_face_colloc(iface)[1];
        double area = grid_.face_area(iface);
        Vector normal = grid_.face_normal(iface);
        double nx = normal.x;
        double ny = normal.y;

        // (u_i + u_j)/2
        {
            double value_x = 0.5 * area * nx;
            double value_y = 0.5 * area * ny;

            if (collocations_.is_boundary_colloc(negative_colloc)) {
                B0.add_value(positive_colloc, negative_colloc, -2 * value_x);
                B1.add_value(positive_colloc, negative_colloc, -2 * value_y);
            } else if (collocations_.is_boundary_colloc(positive_colloc)) {
                B0.add_value(negative_colloc, positive_colloc, 2 * value_x);
                B1.add_value(negative_colloc, positive_colloc, 2 * value_y);
            } else {
                B0.add_value(negative_colloc, negative_colloc, value_x);
                B0.add_value(negative_colloc, positive_colloc, value_x);
                B0.add_value(positive_colloc, negative_colloc, -value_x);
                B0.add_value(positive_colloc, positive_colloc, -value_x);

                B1.add_value(negative_colloc, negative_colloc, value_y);
                B1.add_value(negative_colloc, positive_colloc, value_y);
                B1.add_value(positive_colloc, negative_colloc, -value_y);
                B1.add_value(positive_colloc, positive_colloc, -value_y);
            }
        }

        if (collocations_.is_internal_colloc(negative_colloc) && collocations_.is_internal_colloc(positive_colloc)) {
            // dij * (dp/dn)_ij
            for (const std::pair<const size_t, double>& iter: dfdn_computer_.linear_combination(iface)) {
                size_t column = iter.first;
                double coef = -dn_face_[iface] * area * iter.second;
                B2.add_value(negative_colloc, column, coef);
                B2.add_value(positive_colloc, column, -coef);
            }

            // 0.5 * ( d_i * (dp/dn)_i + d_j * (dp/dn)_ j)
            double d_dpdn_i = dot_product(Vector{grad_p_[negative_colloc].x * du_[negative_colloc],
                                                 grad_p_[negative_colloc].y * dv_[negative_colloc]},
                                          normal);
            double d_dpdn_j = dot_product(Vector{grad_p_[positive_colloc].x * du_[positive_colloc],
                                                 grad_p_[positive_colloc].y * dv_[positive_colloc]},
                                          normal);
            double dpdn = 0.5 * (d_dpdn_i + d_dpdn_j);
            rhs[negative_colloc] += -area * dpdn;
            rhs[positive_colloc] += area * dpdn;
        }
    }

    // boundary condition: dp/dn = 0
    for (size_t icolloc: boundary_info_.bnd_colloc) {
        B0.remove_row(icolloc);
        B1.remove_row(icolloc);
        B2.remove_row(icolloc);
        size_t iface = collocations_.face_index(icolloc);
        double area = grid_.face_area(iface);
        double sgn = (collocations_.tab_face_colloc(iface)[1] == icolloc) ? 1.0 : -1.0;

        for (const std::pair<const size_t, double>& iter: dfdn_computer_.linear_combination(iface)) {
            size_t column = iter.first;
            double coef = sgn * area * iter.second;
            B2.add_value(icolloc, column, coef);
        }
        rhs[icolloc] = 0;
    };
    // p=0 at icell=0
    B0.remove_row(0);
    B1.remove_row(0);
    B2.set_unit_row(0);
    rhs[0] = 0;

    // block assemble
    B = assemble_block_matrix(vec_size(), vec_size(), {{&B0, &B1, &B2}});
};

void CavityFvmCoupledWorker::assemble_slae() {
    std::vector<double> rhs_u(vec_size(), 0);
    std::vector<double> rhs_v(vec_size(), 0);
    std::vector<double> rhs_continuity(vec_size(), 0);
    CsrMatrix Bu, Bv, Bcontinuity;

    u_equation(Bu, rhs_u);
    v_equation(Bv, rhs_v);
    compute_d(Bu, Bv);
    continuity_equation(Bcontinuity, rhs_continuity);

    mat_ = assemble_block_matrix(vec_size(), 3 * vec_size(), {{&Bu}, {&Bv}, {&Bcontinuity}});

    std::copy(rhs_u.begin(), rhs_u.end(), rhs_.begin());
    std::copy(rhs_v.begin(), rhs_v.end(), rhs_.begin() + vec_size());
    std::copy(rhs_continuity.begin(), rhs_continuity.end(), rhs_.begin() + 2 * vec_size());
}

double CavityFvmCoupledWorker::step() {
    AmgcMatrixSolver slv({
        {"solver.maxiter", "1000"},
        {"solver.type", "fgmres"},
        {"precond.coarsening.type", "smoothed_aggregation"},
        {"precond.relax.type", "ilu0"},
    });
    slv.set_matrix(mat_);
    std::vector<double> uvp_old(uvp_);
    slv.solve(rhs_, uvp_);

    // d/dt
    double dudt = 0, dvdt = 0, dpdt = 0;
    for (size_t i = 0; i < vec_size(); ++i) {
        dudt = std::max(dudt, std::abs(uvp_old[i] - uvp_[i]));
        dvdt = std::max(dvdt, std::abs(uvp_old[i + vec_size()] - uvp_[i + vec_size()]));
        dpdt = std::max(dvdt, std::abs(uvp_old[i + 2 * vec_size()] - uvp_[i + 2 * vec_size()]));
    }
    dudt /= delta_t_;
    dvdt /= delta_t_;
    dpdt /= delta_t_;

    // additional values
    grad_p_ = grad_computer_.compute(p_);
    un_face_ = compute_un_rhie_chow();

    to_next_iteration();
    return std::max(std::max(dudt, dvdt), dpdt);
}

void CavityFvmCoupledWorker::save_current_fields(size_t iter) {
    if (writer_) {
        std::string filepath = writer_->add_iter(iter);
        grid_.save_vtk(filepath);
        const size_t nc = grid_.n_cells();
        VtkUtils::add_cell_data(std::vector<double>(p_, p_ + nc), "pressure", filepath, grid_.n_cells());
        VtkUtils::add_cell_vector(std::vector<double>(u_, u_ + nc), std::vector<double>(v_, v_ + nc), "velocity",
                                  filepath, grid_.n_cells());
    }
}

std::vector<double> CavityFvmCoupledWorker::compute_un_rhie_chow() const {

    std::vector<double> ret(grid_.n_faces());
    std::vector<double> dpdn_face = dfdn_computer_.compute(p_);

    for (size_t iface = 0; iface < grid_.n_faces(); ++iface) {
        size_t ci = grid_.tab_face_cell(iface)[0];
        size_t cj = grid_.tab_face_cell(iface)[1];
        if (ci == INVALID_INDEX || cj == INVALID_INDEX) {
            // bc: un = 0 at each boundary
            ret[iface] = 0;
        } else {
            // Rhie-Chow interpolation
            Vector normal = grid_.face_normal(iface);
            Vector uvec_i = Vector(u_[ci], v_[ci]);
            Vector uvec_j = Vector(u_[cj], v_[cj]);

            double u_i = dot_product(uvec_i, normal);
            double u_j = dot_product(uvec_j, normal);

            double d_dpdn_i = dot_product(Vector{grad_p_[ci].x * du_[ci], grad_p_[ci].y * dv_[ci]}, normal);
            double d_dpdn_j = dot_product(Vector{grad_p_[cj].x * du_[cj], grad_p_[cj].y * dv_[cj]}, normal);
            double d_dpdn_ij = dn_face_[iface] * dpdn_face[iface];

            ret[iface] = 0.5 * (u_i + u_j) + 0.5 * (d_dpdn_i + d_dpdn_j - 2 * d_dpdn_ij);
        }
    }

    return ret;
}

TEST_CASE("Cavity FVM-SIMPLE coupled algorithm", "[cavity-fvm-coupled]") {
    std::cout << std::endl << "--- cfd26_test [cavity-fvm-coupled] --- " << std::endl;

    // problem parameters
    double Re = 100;
    size_t max_it = 10'000;
    double eps = 1e-3;
    double tau = 0.5;

    // worker initialization
    RegularGrid2D grid(0, 1, 0, 1, 30, 30);
    // std::string fn = test_directory_file("tetragrid_500.vtk");
    // UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(fn);
    CavityFvmCoupledWorker worker(grid, Re, tau);
    worker.initialize_saver("cavity-fvm-coupled");

    // iterations loop
    size_t it = 0;
    for (it = 1; it < max_it; ++it) {
        double nrm = worker.step();

        // print norm and friction force
        std::cout << it << " " << nrm << std::endl;

        // export solution to vtk
        worker.save_current_fields(it);

        // break if residual is low enough
        if (nrm < eps) {
            break;
        }
    }
    CHECK(it == 20);
}
