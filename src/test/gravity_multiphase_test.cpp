#include "cfd/debug/saver.hpp"
#include "cfd/fvm/fvm_dfdn.hpp"
#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/fvm/fvm_gradient.hpp"
#include "cfd/grid/regular_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd26_test.hpp"

using namespace cfd;

namespace {

const double GX = 0.0;
const double GY = -1.0;
const double PHI_EPS = 1e-6;
const double ALPHA_P = 0.3;
const double ALPHA_U = 0.8;

struct Physics {
    const double rho;
    const double mu;
};

struct Uvphi {
    Uvphi(size_t size, size_t size_faces)
        : u(size, 0),
          v(size, 0),
          phi(size, 0),
          un(size_faces, 0),
          phi_face(size_faces, 0.0) {}

    size_t size() const {
        return u.size();
    }

    // collocated data
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> phi;
    // face data
    std::vector<double> un;
    std::vector<double> phi_face;
};

struct UvSlae {
    double norm_res(const std::vector<double>& u, const std::vector<double>& v, const IGrid& grid) const {
        double ret = 0;
        for (size_t i = 0; i < grid.n_cells(); i++) {
            double res_x = rhs_x[i] - A.mult_vec(i, u);
            double res_y = rhs_y[i] - A.mult_vec(i, v);
            double m = std::max(std::abs(res_x), std::abs(res_y)) / grid.cell_volume(i);
            ret = std::max(ret, m);
        }
        return ret;
    }
    void solve(std::vector<double>& u, std::vector<double>& v) const {
        AmgcMatrixSolver slv;
        slv.set_matrix(A);
        slv.solve(rhs_x, u);
        slv.solve(rhs_y, v);
    }
    size_t size() const {
        return rhs_x.size();
    }
    CsrMatrix A;
    std::vector<double> rhs_x;
    std::vector<double> rhs_y;

    std::vector<double> d; // phi_k * m_i/a_{ii}
    std::vector<double> d_face;
};

class Worker {
public:
    Worker(Physics physics1, Physics physics2, std::shared_ptr<IGrid> grid, double dt)
        : phys1_(physics1),
          phys2_(physics2),
          dt_(dt),
          grid_(grid),
          ecolloc_(*grid_),
          grad_computer_(*grid_, ecolloc_),
          dfdn_computer_(*grid_, ecolloc_),
          sol1_(ecolloc_.size(), grid_->n_faces()),
          sol2_(ecolloc_.size(), grid_->n_faces()),
          sol1_old_(ecolloc_.size(), grid_->n_faces()),
          sol2_old_(ecolloc_.size(), grid_->n_faces()) {

        p_.resize(ecolloc_.size(), 0);
        // Initial condition
        for (size_t i = 0; i < ecolloc_.points.size(); i++) {
            if (ecolloc_.points[i].y <= 0.5) {
                sol1_.phi[i] = PHI_EPS;
            } else {
                sol1_.phi[i] = 1.0 - PHI_EPS;
            }
            sol2_.phi[i] = 1.0 - sol1_.phi[i];
        }
    }

    void init_time_step() {
        sol1_old_ = sol1_;
        sol2_old_ = sol2_;
        // Compute phi1, phi2
        sol1_.phi = compute_phi_explicit(sol1_old_.un, sol1_old_.phi);
        for (size_t i = 0; i < sol1_.size(); i++) {
            sol1_.phi[i] = std::clamp(sol1_.phi[i], 0.0 + PHI_EPS, 1.0 - PHI_EPS);
            sol2_.phi[i] = 1 - sol1_.phi[i];
        }
        for (size_t i = 0; i < grid_->n_faces(); ++i) {
            sol1_.phi_face[i] = ecolloc_.face_approx(i, sol1_.phi);
            sol2_.phi_face[i] = ecolloc_.face_approx(i, sol2_.phi);
        }
    }

    double init_simple_step() {
        // Momentum matrix assemble
        slae1_ = assemble_momentum_slae(sol1_, sol1_old_, phys1_);
        slae2_ = assemble_momentum_slae(sol2_, sol2_old_, phys2_);

        // compute norms
        double r1 = slae1_.norm_res(sol1_.u, sol1_.v, *grid_);
        double r2 = slae2_.norm_res(sol2_.u, sol2_.v, *grid_);

        return std::max(r1, r2);
    }

    void step() {
        // 1. Ustar
        std::vector<double> u1star(sol1_.u), u2star(sol1_.u), v1star(sol1_.u), v2star(sol1_.u);
        slae1_.solve(u1star, v1star);
        slae2_.solve(u2star, v2star);
        // 2. Un_star
        std::vector<double> un1_face_star = compute_un_face_rhie_chow(slae1_, u1star, v1star);
        std::vector<double> un2_face_star = compute_un_face_rhie_chow(slae2_, u2star, v2star);
        // 3. p_prime
        std::vector<double> p_prime = compute_pressure(un1_face_star, un2_face_star);
        std::vector<Vector> grad_p_prime = grad_computer_.compute(p_prime);
        // 4. U'
        auto [u1_prime, v1_prime] = compute_u_prime(slae1_, grad_p_prime);
        auto [u2_prime, v2_prime] = compute_u_prime(slae2_, grad_p_prime);
        // 5. Assemble current values
        for (size_t i = 0; i < ecolloc_.size(); ++i) {
            p_[i] += ALPHA_P * p_prime[i];
            sol1_.u[i] = u1star[i] + u1_prime[i];
            sol1_.v[i] = v1star[i] + v1_prime[i];
            sol2_.u[i] = u2star[i] + u2_prime[i];
            sol2_.v[i] = v2star[i] + v2_prime[i];
        }
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            double dp_prime_dn = dfdn_computer_.compute(iface, p_prime);
            double un1_face_prime = -slae1_.d_face[iface] * dp_prime_dn;
            double un2_face_prime = -slae2_.d_face[iface] * dp_prime_dn;

            sol1_.un[iface] = un1_face_star[iface] + un1_face_prime;
            sol2_.un[iface] = un2_face_star[iface] + un2_face_prime;
        }
    }

    void save_solution(VtkUtils::TimeSeriesWriter& writer, double time) const {
        std::string fname = writer.add(time);
        if (!fname.empty()) {
            grid_->save_vtk(fname);
            VtkUtils::add_cell_data(sol1_.phi, "phi1", fname, grid_->n_cells());
            VtkUtils::add_cell_data(sol2_.phi, "phi2", fname, grid_->n_cells());
            VtkUtils::add_cell_data(p_, "p", fname, grid_->n_cells());

            VtkUtils::add_cell_vector(sol1_.u, sol1_.v, "u1", fname, grid_->n_cells());
            VtkUtils::add_cell_vector(sol2_.u, sol2_.v, "u2", fname, grid_->n_cells());
        }
    }

private:
    const Physics phys1_;
    const Physics phys2_;
    const double dt_;
    const std::shared_ptr<IGrid> grid_;
    const FvmExtendedCollocations ecolloc_;
    const GaussLinearFvmCellGradient grad_computer_;
    const FvmFacesDn dfdn_computer_;

    Uvphi sol1_;
    Uvphi sol2_;
    std::vector<double> p_;
    Uvphi sol1_old_;
    Uvphi sol2_old_;

    UvSlae slae1_;
    UvSlae slae2_;

    std::vector<double> compute_phi_explicit(const std::vector<double>& un, const std::vector<double>& phi_old) const {
        std::vector<double> phi = phi_old;

        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            double flux = un[iface] * grid_->face_area(iface);
            auto [n_colloc, p_colloc] = ecolloc_.tab_face_colloc(iface);
            if (flux < 0) {
                std::swap(n_colloc, p_colloc);
                flux *= -1;
            }
            // compute upwind value
            double face_value = phi_old[n_colloc] * flux * dt_;
            // add to values
            if (ecolloc_.is_internal_colloc(n_colloc)) {
                phi[n_colloc] -= face_value / grid_->cell_volume(n_colloc);
            }
            if (ecolloc_.is_internal_colloc(p_colloc)) {
                phi[p_colloc] += face_value / grid_->cell_volume(p_colloc);
            }
        }

        return phi;
    }

    UvSlae assemble_momentum_slae(const Uvphi& sol, const Uvphi& sol_old, const Physics& phys) const {
        std::vector<double> rhs_x(ecolloc_.size(), 0);
        std::vector<double> rhs_y(ecolloc_.size(), 0);

        LodMatrix mat(ecolloc_.size());
        // =============== Cell loop
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            double volume = grid_->cell_volume(icell);
            // rho * d (phi * u) / dt
            double vl = phys.rho * volume / dt_ * sol.phi[icell];
            double vr = phys.rho * volume / dt_ * sol_old.phi[icell];
            mat.add_value(icell, icell, vl);
            rhs_x[icell] += vr * sol_old.u[icell];
            rhs_y[icell] += vr * sol_old.v[icell];

            // = -phi * (grad p)
            Vector grad_p = grad_computer_.cell_compute(icell, p_);
            rhs_x[icell] -= sol.phi[icell] * volume * grad_p.x;
            rhs_y[icell] -= sol.phi[icell] * volume * grad_p.y;

            // phi * rho * g
            rhs_x[icell] += volume * sol.phi[icell] * phys.rho * GX;
            rhs_y[icell] += volume * sol.phi[icell] * phys.rho * GY;
        }

        // =============== Faces loop
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [n_colloc, p_colloc] = ecolloc_.tab_face_colloc(iface);
            double area = grid_->face_area(iface);
            double phi_face = sol.phi_face[iface];

            // Convection: div( rho u @ u)
            {
                // Symmetric
                double coef = phys.rho * area * sol.un[iface] * phi_face / 2.0;
                mat.add_value(n_colloc, n_colloc, coef);
                mat.add_value(n_colloc, p_colloc, coef);
                mat.add_value(p_colloc, p_colloc, -coef);
                mat.add_value(p_colloc, n_colloc, -coef);
            }

            // Diffusion: - div(phi*mu*grad(u))  including (-mu du/dn) approx on boundary collocations
            for (const auto& [column, coef]: dfdn_computer_.linear_combination(iface)) {
                double flux = -phys.mu * area * coef * phi_face;
                mat.add_value(n_colloc, column, flux);
                mat.add_value(p_colloc, column, -flux);
            }
        }
        // relax
        for (size_t irow = 0; irow < grid_->n_cells(); ++irow) {
            double diag = mat.value(irow, irow);
            mat.set_value(irow, irow, diag / ALPHA_U);
            rhs_x[irow] += diag / ALPHA_U * (1 - ALPHA_U) * sol.u[irow];
            rhs_y[irow] += diag / ALPHA_U * (1 - ALPHA_U) * sol.v[irow];
        }

        // compute d = (phi * m_i / a_ii)
        std::vector<double> d(ecolloc_.size(), 0);
        for (size_t i = 0; i < grid_->n_cells(); i++) {
            d[i] = sol.phi[i] * grid_->cell_volume(i) / mat.value(i, i);
        }
        for (size_t i = grid_->n_cells(); i < d.size(); i++) {
            size_t o = ecolloc_.face_index(i);
            auto [n, p] = grid_->tab_face_cell(o);
            d[i] = (n != INVALID_INDEX) ? d[n] : d[p];
        }
        std::vector<double> dface(grid_->n_faces(), 0);
        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            dface[iface] = ecolloc_.face_approx(iface, d);
        }
        // boundary conditions: u=0 on all b.faces
        for (size_t icolloc: ecolloc_.face_collocations) {
            mat.set_unit_row(icolloc);
            rhs_x[icolloc] = 0;
            rhs_y[icolloc] = 0;
        }

        return UvSlae{mat.to_csr(), std::move(rhs_x), std::move(rhs_y), std::move(d), std::move(dface)};
    }

    std::vector<double> compute_un_face_rhie_chow(const UvSlae& slae, const std::vector<double>& u_star,
                                                  const std::vector<double>& v_star) const {
        std::vector<double> ret(grid_->n_faces());
        std::vector<Vector> grad_p = grad_computer_.compute(p_);

        for (size_t iface = 0; iface < grid_->n_faces(); ++iface) {
            auto [ci, cj] = grid_->tab_face_cell(iface);
            double dpdn_face = dfdn_computer_.compute(iface, p_);
            if (ci == INVALID_INDEX || cj == INVALID_INDEX) {
                // boundary condition: u_n = 0
                ret[iface] = 0;
            } else {
                // Rhie-Chow interpolation
                Vector normal = grid_->face_normal(iface);
                Vector uvec_i = Vector(u_star[ci], v_star[ci]);
                Vector uvec_j = Vector(u_star[cj], v_star[cj]);

                double ustar_i = dot_product(uvec_i, normal);
                double ustar_j = dot_product(uvec_j, normal);
                double d_dpdn_i = dot_product(Vector{grad_p[ci].x * slae.d[ci], grad_p[ci].y * slae.d[ci]}, normal);
                double d_dpdn_j = dot_product(Vector{grad_p[cj].x * slae.d[cj], grad_p[cj].y * slae.d[cj]}, normal);
                double d_dpdn_ij = slae.d_face[iface] * dpdn_face;

                ret[iface] = 0.5 * (ustar_i + ustar_j) + 0.5 * (d_dpdn_i + d_dpdn_j - 2 * d_dpdn_ij);
            }
        }

        return ret;
    }
    std::vector<double> compute_pressure(const std::vector<double>& un_face_star1,
                                         const std::vector<double>& un_face_star2) const {
        LodMatrix mat(ecolloc_.size());
        std::vector<double> rhs(ecolloc_.size(), 0);

        for (size_t iface = 0; iface < grid_->n_faces(); iface++) {
            double area = grid_->face_area(iface);
            auto [i, j] = ecolloc_.tab_face_colloc(iface);
            // diffusion (left)
            double d1 = slae1_.d_face[iface];
            double d2 = slae2_.d_face[iface];
            double phi1 = sol1_.phi_face[iface];
            double phi2 = sol2_.phi_face[iface];
            double d_ij = phi1 * d1 + phi2 * d2;
            for (auto [icol, coef]: dfdn_computer_.linear_combination(iface)) {
                mat.add_value(i, icol, -coef * d_ij * area);
                mat.add_value(j, icol, coef * d_ij * area);
            }
            double un_star =
                un_face_star1[iface] * sol1_.phi_face[iface] + un_face_star2[iface] * sol2_.phi_face[iface];
            // right
            if (!ecolloc_.is_boundary_colloc(i)) {
                rhs[i] -= un_star * area;
            }
            if (!ecolloc_.is_boundary_colloc(j)) {
                rhs[j] += un_star * area;
            }
        }
        mat.set_unit_row(0);
        rhs[0] = 0;
        std::vector<double> pressure;
        AmgcMatrixSolver::solve_slae(mat.to_csr(), rhs, pressure);
        return pressure;
    }

    std::pair<std::vector<double>, std::vector<double>> compute_u_prime(const UvSlae& slae,
                                                                        const std::vector<Vector>& grad_p) const {
        std::vector<double> u(slae.size(), 0);
        std::vector<double> v(slae.size(), 0);

        for (size_t i = 0; i < grid_->n_cells(); ++i) {
            u[i] = -slae.d[i] * grad_p[i].x;
            v[i] = -slae.d[i] * grad_p[i].y;
        }

        return {u, v};
    }
};
} // namespace

TEST_CASE("Multiphase gravity segragation", "[gravity-mf]") {
    std::cout << std::endl << "--- cfd_test [gravity-mf] --- " << std::endl;

    Physics water{1.0, 1.0 / 200.0};
    Physics oil{0.8, 1.0 / 100.0};

    double dt = 2e-2;
    double end_time = 0.1;
    double simple_eps = 1e-3;
    size_t max_simple_steps = 50;

    auto grid = std::make_shared<RegularGrid2D>(0, 1, 0, 1, 30, 30);

    Worker worker(water, oil, grid, dt);

    VtkUtils::TimeSeriesWriter writer("gravity-mf");
    writer.set_time_step(0.5);

    double time = 0;
    double nrm;
    size_t it;
    worker.save_solution(writer, 0);
    while (time + dt < end_time + 1e-6) {
        time += dt;
        std::cout << "time=" << time << std::endl;
        worker.init_time_step();
        for (it = 0; true; ++it) {
            nrm = worker.init_simple_step();
            std::cout << "--iter: " << it << "; residual = " << nrm << std::endl;
            worker.step();
            if (it > 0) {
                if (nrm < simple_eps) {
                    std::cout << "converged in " << it << " steps" << std::endl;
                    break;
                } else if (it >= max_simple_steps) {
                    std::cout << "WARNING: failed to converge with residual = " << nrm << std::endl;
                    break;
                }
            }
        }
        worker.save_solution(writer, time);
    }
    CHECK(nrm == Approx(0.0007046409));
    CHECK(it == 7);
}
