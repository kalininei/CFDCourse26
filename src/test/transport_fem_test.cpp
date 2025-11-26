#include "cfd/fem/elem1d/segment_linear.hpp"
#include "cfd/fem/elem2d/triangle_linear.hpp"
#include "cfd/fem/fem_assembler.hpp"
#include "cfd/grid/cell_finder.hpp"
#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/matrix_iter.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd/numeric_integration/numeric_integration.hpp"
#include "cfd26_test.hpp"
#include "test/utils/filesystem.hpp"

using namespace cfd;

namespace {

////////////////////////////////////////////
// Solution
////////////////////////////////////////////
class ISolution {
public:
    virtual ~ISolution() = default;

    virtual Vector velocity(Point p) const = 0;
    virtual double init_solution(Point p) const = 0;
    virtual double exact_solution(Point p, double t) const = 0;
};

class Solution1D : public ISolution {
public:
    Vector velocity(Point) const override {
        return Vector{1};
    }
    double init_solution(Point p) const override {
        constexpr double sigma = 0.1;
        return exp(-p.x * p.x / sigma / sigma);
    }
    double exact_solution(Point p, double t) const override {
        return init_solution(p - t * Vector{1});
    }
};

//////////////////////////////////////
// IFemBuilder
//////////////////////////////////////
class IFemBuilder {
public:
    virtual ~IFemBuilder() = default;

    virtual FemAssembler build() const = 0;
    virtual std::vector<size_t> boundary_bases() const = 0;
};

//////////////////////////////////////
// FemLinearSegment
//////////////////////////////////////
class FemLinearSegment : public IFemBuilder {
public:
    FemLinearSegment(std::shared_ptr<Grid1D> grid) : grid_(grid) {}

    FemAssembler build() const override {
        size_t n_bases = grid_->n_points();
        std::vector<FemElement> elements;
        std::vector<std::vector<size_t>> tab_elem_basis;

        // elements
        for (size_t icell = 0; icell < grid_->n_cells(); ++icell) {
            auto ipoints = grid_->tab_cell_point(icell);
            if (ipoints.size() != 2) {
                // only segment grid is allowed
                _THROW_NOT_IMP_;
            }
            Point p0 = grid_->point(ipoints[0]);
            Point p1 = grid_->point(ipoints[1]);

            auto geom = std::make_shared<SegmentLinearGeometry>(p0, p1);
            auto basis = std::make_shared<SegmentLinearBasis>();
            auto quad = quadrature_segment_gauss3();
            FemElement elem{geom, basis, quad};

            elements.push_back(elem);
            std::vector<size_t> tab = ipoints;
            tab_elem_basis.push_back(tab);
        }

        return FemAssembler(n_bases, elements, tab_elem_basis);
    }

    std::vector<size_t> boundary_bases() const override {
        return grid_->boundary_points();
    }

private:
    std::shared_ptr<Grid1D> grid_;
};
////////////////////////////////////////////////////////////
// Basic Worker
////////////////////////////////////////////////////////////

class ATestTransportWorker {
public:
    virtual ~ATestTransportWorker() {}

    Vector velocity(Point p) const {
        return solution_.velocity(p);
    }
    double init_solution(Point p) const {
        return solution_.init_solution(p);
    }
    double exact_solution(Point p) const {
        return solution_.exact_solution(p, time_);
    }

    ATestTransportWorker(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : grid_(grid), fem_(builder.build()), tau_(tau), u_(fem_.n_bases(), 0.0), solution_(solution),
          boundary_bases_(builder.boundary_bases()) {

        // mass matrix
        mass_.set_stencil(fem_.stencil());
        for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
            std::vector<double> local_mass = element_mass_matrix(ielem);
            fem_.add_to_global_matrix(ielem, local_mass, mass_.vals());
        }
        // transport matrix
        high_order_transport_.set_stencil(fem_.stencil());
        for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
            std::vector<double> local_transport = element_transport_matrix(ielem);
            fem_.add_to_global_matrix(ielem, local_transport, high_order_transport_.vals());
        }

        // load vector
        load_vector_.resize(fem_.n_bases(), 0);
        for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
            std::vector<double> local_load = element_load_vector(ielem);
            fem_.add_to_global_vector(ielem, local_load, load_vector_);
        }

        // diffusion: dij = max(0, kij, kji), dii = sum_j (dij)
        transport_diffusion_.set_stencil(fem_.stencil());
        for (size_t irow = 0; irow < fem_.n_bases(); ++irow) {
            size_t a_start = fem_.stencil().addr()[irow];
            size_t a_end = fem_.stencil().addr()[irow + 1];
            double diag = 0.0;
            size_t a_diag = INVALID_INDEX;
            for (size_t a = a_start; a < a_end; ++a) {
                size_t jcol = fem_.stencil().cols()[a];
                if (jcol != irow) {
                    double kij = high_order_transport_.vals()[a];
                    double kji = high_order_transport_.value(jcol, irow);
                    double dij = std::max(std::max(0.0, kij), kji);
                    transport_diffusion_.vals()[a] = dij;
                    diag -= dij;
                } else {
                    a_diag = a;
                }
            }
            if (a_diag == INVALID_INDEX) {
                _THROW_INTERNAL_ERROR_;
            }
            transport_diffusion_.vals()[a_diag] = diag;
        }

        // _low order transport: L = K + D
        low_order_transport_.set_stencil(fem_.stencil());
        for (size_t i = 0; i < fem_.stencil().n_nonzeros(); ++i) {
            low_order_transport_.vals()[i] = high_order_transport_.vals()[i] + transport_diffusion_.vals()[i];
        }

        // initial solution
        for (size_t ibas = 0; ibas < u_.size(); ++ibas) {
            Point p = fem_.reference_point(ibas);
            u_[ibas] = init_solution(p);
        }
    }

    void step() {
        time_ += tau_;
        impl_step();
    }

    void save_vtk(const std::string& filename) {
        // save grid
        grid_->save_vtk(filename);

        // save numerical solution
        VtkUtils::add_point_data(u_, "numerical", filename, grid_->n_points());

        // save exact solution
        std::vector<double> exact(grid_->n_points());
        for (size_t i = 0; i < grid_->n_points(); ++i) {
            exact[i] = exact_solution(grid_->point(i));
        }
        VtkUtils::add_point_data(exact, "exact", filename);
    }

    double current_time() const {
        return time_;
    }

    std::vector<Vector> velocity() const {
        std::vector<Vector> ret;
        for (size_t ibas = 0; ibas < fem_.n_bases(); ++ibas) {
            Point p = fem_.reference_point(ibas);
            ret.push_back(solution_.velocity(p));
        }
        return ret;
    }

    double compute_norm2() {
        std::vector<double> diff(u_.size());
        for (size_t ibas = 0; ibas < fem_.n_bases(); ++ibas) {
            Point p = fem_.reference_point(ibas);
            diff[ibas] = u_[ibas] - exact_solution(p);
        }
        return compute_norm2(diff);
    }

    double compute_norm2(const std::vector<double>& v) {
        double integral = 0;
        double full_area = 0;
        for (size_t ibas = 0; ibas < fem_.n_bases(); ++ibas) {
            integral += load_vector_[ibas] * (v[ibas] * v[ibas]);
            full_area += load_vector_[ibas];
        }
        return std::sqrt(integral / full_area);
    }

    const CsrMatrix& high_order_transport() const {
        return high_order_transport_;
    }
    const CsrMatrix& low_order_transport() const {
        return low_order_transport_;
    }
    const CsrMatrix& transport_diffusion() const {
        return transport_diffusion_;
    }
    const std::vector<double>& load_vector() const {
        return load_vector_;
    }
    const FemAssembler& fem() const {
        return fem_;
    }
    const IGrid& grid() const {
        return *grid_;
    }

protected:
    std::shared_ptr<IGrid> grid_;
    FemAssembler fem_;
    CsrMatrix mass_;
    std::vector<double> load_vector_;
    CsrMatrix high_order_transport_;
    CsrMatrix transport_diffusion_;
    CsrMatrix low_order_transport_;

    const double tau_;
    std::vector<double> u_;
    double time_ = 0;
    const ISolution& solution_;
    std::vector<size_t> boundary_bases_;

    const std::vector<size_t>& boundary_bases() const {
        return boundary_bases_;
    }

    void apply_dirichlet_bc(std::vector<double>& s) const {
        for (size_t ibas: boundary_bases()) {
            Point p = fem_.reference_point(ibas);
            s[ibas] = exact_solution(p);
        }
    }
    void apply_uniform_dirichlet_bc(std::vector<double>& s) const {
        for (size_t ibas: boundary_bases()) {
            s[ibas] = 0;
        }
    }
    void apply_dirichlet_bc(CsrMatrix& m) const {
        for (size_t ibas: boundary_bases()) {
            m.set_unit_row(ibas);
        }
    }

private:
    virtual void impl_step() = 0;

    std::vector<double> element_load_vector(size_t ielem) const {
        const FemElement& el = fem_.element(ielem);

        auto fun = [&el](Point p) -> std::vector<double> {
            std::vector<double> ret(el.basis->size(), 0.0);
            auto val = FemElementValue(&el, p);
            for (size_t ibas = 0; ibas < ret.size(); ++ibas) {
                ret[ibas] = val.phi(ibas) * val.modj();
            }
            return ret;
        };

        return el.quadrature->integrate(fun);
    }

    std::vector<double> element_mass_matrix(size_t ielem) const {
        const FemElement& el = fem_.element(ielem);

        auto fun = [&el](Point p) -> std::vector<double> {
            const size_t n = el.basis->size();
            std::vector<double> ret(n * n, 0.0);
            auto val = FemElementValue(&el, p);

            for (size_t ibas = 0; ibas < n; ++ibas) {
                for (size_t jbas = 0; jbas <= ibas; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    ret[k1] = val.phi(ibas) * val.phi(jbas) * val.modj();
                    if (ibas != jbas) {
                        const size_t k2 = jbas * n + ibas;
                        ret[k2] = ret[k1];
                    }
                }
            }

            return ret;
        };

        return el.quadrature->integrate(fun);
    }

    std::vector<double> element_transport_matrix(size_t ielem) const {
        const FemElement& el = fem_.element(ielem);
        std::vector<Vector> v1 = velocity();

        auto fun = [&](Point p) -> std::vector<double> {
            const size_t n = el.basis->size();
            const std::vector<Vector> element_velocity = fem_.local_vector(ielem, v1);

            std::vector<double> ret(n * n, 0.0);
            auto val = FemElementValue(&el, p);

            for (size_t ibas = 0; ibas < n; ++ibas) {
                for (size_t jbas = 0; jbas < n; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    Vector vel = val.interpolate(element_velocity);
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] = -dot_product(vel, g2) * val.phi(ibas) * val.modj();
                }
            }

            return ret;
        };

        return el.quadrature->integrate(fun);
    }
};

} // namespace

///////////////////////////////////////////////////////////////////////////////
// Linear transport
///////////////////////////////////////////////////////////////////////////////

namespace {

class ExplicitUpwind : public ATestTransportWorker {
public:
    ExplicitUpwind(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : ATestTransportWorker(grid, tau, builder, solution) {}

private:
    void impl_step() override {
        std::vector<double> u_new(u_.size(), 0);
        for (size_t i = 0; i < fem_.n_bases(); ++i) {
            double h = load_vector_[i];
            double lu = low_order_transport_.mult_vec(i, u_);
            u_new[i] = u_[i] + tau_ / h * lu;
        }
        apply_dirichlet_bc(u_new);
        std::swap(u_, u_new);
    }
};

class CrankNicolsonUpwind : public ATestTransportWorker {
public:
    CrankNicolsonUpwind(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : ATestTransportWorker(grid, tau, builder, solution), solver_(1000, 1e-12) {
        // LHS = MASS - 0.5 * tau * L
        CsrMatrix lhs(fem_.stencil());
        for (size_t i = 0; i < lhs.n_nonzeros(); ++i) {
            lhs.vals()[i] = mass_.vals()[i] - 0.5 * tau_ * low_order_transport_.vals()[i];
        }
        // Dirichlet BC
        apply_dirichlet_bc(lhs);

        solver_.set_matrix(lhs);
    }

private:
    AmgcMatrixSolver solver_;

    void impl_step() override {
        // RHS = (MASS + 0.5 * tau * L) * u_old
        std::vector<double> rhs(fem_.n_bases(), 0);
        for (size_t i = 0; i < fem_.n_bases(); ++i) {
            rhs[i] = mass_.mult_vec(i, u_) + tau_ * 0.5 * low_order_transport_.mult_vec(i, u_);
        }
        // Dirichlet bc
        apply_dirichlet_bc(rhs);

        // solver
        solver_.solve(rhs, u_);
    }
};

class CrankNicolsonCentral : public ATestTransportWorker {
public:
    CrankNicolsonCentral(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : ATestTransportWorker(grid, tau, builder, solution), solver_(1000, 1e-12) {
        // LHS = MASS - 0.5 * tau * K
        CsrMatrix lhs(fem_.stencil());
        for (size_t i = 0; i < lhs.n_nonzeros(); ++i) {
            lhs.vals()[i] = mass_.vals()[i] - 0.5 * tau_ * high_order_transport_.vals()[i];
        }
        // Dirichlet BC
        apply_dirichlet_bc(lhs);

        solver_.set_matrix(lhs);
    }

private:
    AmgcMatrixSolver solver_;

    void impl_step() override {
        // RHS = (MASS + 0.5 * tau * K) * u_old
        std::vector<double> rhs(fem_.n_bases(), 0);
        for (size_t i = 0; i < fem_.n_bases(); ++i) {
            rhs[i] = mass_.mult_vec(i, u_) + tau_ * 0.5 * high_order_transport_.mult_vec(i, u_);
        }
        // Dirichlet bc
        apply_dirichlet_bc(rhs);

        // solver
        solver_.solve(rhs, u_);
    }
};

class AFct : public ATestTransportWorker {
public:
    AFct(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : ATestTransportWorker(grid, tau, builder, solution) {
        mass_lumped_ = load_vector_;
    }

protected:
    std::vector<double> mass_lumped_;

    std::vector<double> compute_low_order_solution(double theta) {
        std::vector<double> ret(u_);
        for (size_t i = 0; i < fem_.n_bases(); ++i) {
            ret[i] += (1 - theta) * tau_ / mass_lumped_[i] * low_order_transport_.mult_vec(i, u_);
        }
        apply_dirichlet_bc(ret);
        return ret;
    }

    std::vector<double> compute_corrected_antidiffusion(const std::vector<double>& low_order_u,
                                                        const std::vector<double>& u_next, double theta) {
        CsrMatrix f = compute_antidiffusion(u_next, theta);

        for (size_t irow = 0; irow < f.n_rows(); ++irow) {
            size_t a0 = f.addr()[irow];
            size_t a1 = f.addr()[irow + 1];
            for (size_t a = a0; a < a1; ++a) {
                size_t jcol = f.cols()[a];
                if (irow != jcol) {
                    double val = f.vals()[a] * (low_order_u[jcol] - low_order_u[irow]);
                    if (val > 0) {
                        f.vals()[a] = 0;
                    }
                }
            }
        }

        use_zalesaks_limiter(low_order_u, f);
        std::vector<double> ret(low_order_u.size(), 0);
        for (size_t i = 0; i < ret.size(); ++i) {
            ret[i] = f.row_sum(i);
        }
        return ret;
    }

    CsrMatrix compute_antidiffusion(const std::vector<double>& u_next, double theta) {
        CsrMatrix ret(fem_.stencil());
        // dudt
        std::vector<double> dudt(u_next.size());
        for (size_t i = 0; i < dudt.size(); ++i) {
            dudt[i] = (u_next[i] - u_[i]) / tau_;
        }
        // fij = mij * (dudt_i - dudt_j)
        //       + theta     * dij * (u_next_i - u_next_j)
        //       + (1-theta) * dij * (u_i - u_j)
        const std::vector<size_t>& addr = ret.addr();
        const std::vector<size_t>& cols = ret.cols();
        std::vector<double>& f = ret.vals();
        std::vector<double>& m = mass_.vals();
        std::vector<double>& d = transport_diffusion_.vals();

        for (size_t irow = 0; irow < ret.n_rows(); ++irow) {
            for (size_t a = addr[irow]; a < addr[irow + 1]; ++a) {
                size_t jcol = cols[a];
                if (irow != jcol) {
                    f[a] = m[a] * (dudt[irow] - dudt[jcol]) + theta * d[a] * (u_next[irow] - u_next[jcol]) +
                           (1 - theta) * d[a] * (u_[irow] - u_[jcol]);
                }
            }
        }

        return ret;
    }

    void use_zalesaks_limiter(const std::vector<double>& low_order_u, CsrMatrix& f) {
        const std::vector<size_t>& addr = f.addr();
        const std::vector<size_t>& cols = f.cols();
        std::vector<double>& fvals = f.vals();

        // compute r +/-
        std::vector<double> r_plus(f.n_rows(), 0.0);
        std::vector<double> r_minus(f.n_rows(), 0.0);
        for (size_t irow = 0; irow < f.n_rows(); ++irow) {
            double p_plus = 0;
            double p_minus = 0;
            double q_plus = 0;
            double q_minus = 0;
            for (size_t a = addr[irow]; a < addr[irow + 1]; ++a) {
                size_t jcol = cols[a];
                if (jcol != irow) {
                    p_plus += std::max(0.0, fvals[a]);
                    p_minus += std::min(0.0, fvals[a]);
                    q_plus = std::max(q_plus, low_order_u[jcol] - low_order_u[irow]);
                    q_minus = std::min(q_minus, low_order_u[jcol] - low_order_u[irow]);
                }
            }
            if (std::abs(p_plus) < 1e-12) {
                p_plus = 1e-12;
            }
            if (std::abs(p_minus) < 1e-12) {
                p_minus = 1e-12;
            }
            r_plus[irow] = std::min(1.0, (mass_lumped_[irow] * q_plus) / (tau_ * p_plus));
            r_minus[irow] = std::min(1.0, (mass_lumped_[irow] * q_minus) / (tau_ * p_minus));
        }

        // use limiter
        for (size_t irow = 0; irow < f.n_rows(); ++irow) {
            for (size_t a = addr[irow]; a < addr[irow + 1]; ++a) {
                size_t jcol = cols[a];
                if (jcol != irow) {
                    double aij =
                        (fvals[a] > 0) ? std::min(r_plus[irow], r_minus[jcol]) : std::min(r_minus[irow], r_plus[jcol]);
                    fvals[a] *= aij;
                }
            }
        }
    }
};

class ExplicitFct : public AFct {
    static constexpr size_t MAX_ITERS = 1000;
    static constexpr double ITERS_EPS = 1e-8;

public:
    ExplicitFct(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : AFct(grid, tau, builder, solution) {}

private:
    void impl_step() override {
        // 1. Low order solution
        auto low_order_u = compute_low_order_solution(0.0);

        std::vector<double> u_new = low_order_u;
        for (size_t it = 0; it < MAX_ITERS; ++it) {
            // 2. compute antidiffusion fluxes
            std::vector<double> antidiffusion = compute_corrected_antidiffusion(low_order_u, u_new, 0.0);

            // 3. compute residual
            std::vector<double> residual(u_.size(), 0);
            for (size_t i = 0; i < u_.size(); ++i) {
                residual[i] = mass_lumped_[i] * low_order_u[i] + tau_ * antidiffusion[i] - mass_lumped_[i] * u_new[i];
            }
            apply_uniform_dirichlet_bc(residual);
            double n2 = compute_norm2(residual);
            if (n2 < ITERS_EPS) {
                break;
            } else if (it == MAX_ITERS - 1) {
                std::cout << "WARNING: FCT loop failed to converge with n2=" << n2 << std::endl;
            }
            for (size_t i = 0; i < u_.size(); ++i) {
                double delta_u = residual[i] / mass_lumped_[i];
                u_new[i] += 0.99 * delta_u; // use stabilized multiplier = 0.99 for better convergence
            }
        }
        u_ = u_new;
    }
};

class CrankNicolsonFct : public AFct {
    static constexpr size_t MAX_ITERS = 1000;
    static constexpr double ITERS_EPS = 1e-8;

public:
    CrankNicolsonFct(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : AFct(grid, tau, builder, solution), solver_(1000, 1e-12) {
        // LHS = MASS_LUMPED - 0.5 * tau * L
        lhs_.set_stencil(fem_.stencil());
        lhs_.set_diagonal(mass_lumped_);
        for (size_t i = 0; i < lhs_.n_nonzeros(); ++i) {
            lhs_.vals()[i] -= 0.5 * tau_ * low_order_transport_.vals()[i];
        }
        // Dirichlet BC
        apply_dirichlet_bc(lhs_);
        solver_.set_matrix(lhs_);
    }

private:
    AmgcMatrixSolver solver_;
    CsrMatrix lhs_;

    void impl_step() override {
        // 1. Low order solution
        auto low_order_u = compute_low_order_solution(0.5);

        std::vector<double> u_new = low_order_u;
        for (size_t it = 0; it < MAX_ITERS; ++it) {
            // 2. compute antidiffusion fluxes
            std::vector<double> antidiffusion = compute_corrected_antidiffusion(low_order_u, u_new, 0.5);

            // 3. compute residual
            std::vector<double> residual(u_.size(), 0);
            for (size_t i = 0; i < u_.size(); ++i) {
                residual[i] = mass_lumped_[i] * low_order_u[i] + tau_ * antidiffusion[i] - lhs_.mult_vec(i, u_new);
            }
            apply_uniform_dirichlet_bc(residual);
            double n2 = compute_norm2(residual);
            if (n2 < ITERS_EPS) {
                break;
            } else if (it == MAX_ITERS - 1) {
                std::cout << "WARNING: FCT loop failed to converge with n2=" << n2 << std::endl;
            }
            std::vector<double> delta_u(u_.size(), 0);
            solver_.solve(residual, delta_u);
            for (size_t i = 0; i < u_.size(); ++i) {
                u_new[i] += delta_u[i];
            }
        }
        u_ = u_new;
    }
};

} // namespace

TEST_CASE("Transport 1D solver, simple fem schemes", "[transport1-fem]") {
    std::cout << std::endl << "--- cfd_test [transport1-fem] --- " << std::endl;
    // parameters
    const double tend = 0.5;
    const double tau = 1e-2;
    const double save_tau = 0.05;
    const size_t n_cells = 20;

    auto compute = [&](ATestTransportWorker& worker, std::string out_file) {
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

    Solution1D solution;
    auto grid = std::make_shared<Grid1D>(0, 1, n_cells);
    FemLinearSegment builder(grid);

    // =========== EXPLICIT worker
    {
        // UPWIND
        ExplicitUpwind worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.193719381).margin(1e-5));
    }

    // ============ Crank-Nicolson
    {
        // UPWIND
        CrankNicolsonUpwind worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.2093664244).margin(1e-5));
    }
    {
        // CENTRAL
        CrankNicolsonCentral worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.0249143615).margin(1e-5));
    }
}

TEST_CASE("Transport 1D solver, fct fem scheme", "[transport1-fem-fct]") {
    std::cout << std::endl << "--- cfd_test [transport1-fem-fct] --- " << std::endl;
    // parameters
    const double tend = 0.5;
    const double tau = 1e-2;
    const double save_tau = 0.05;
    const size_t n_cells = 20;

    auto compute = [&](ATestTransportWorker& worker, std::string out_file) {
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

    Solution1D solution;
    auto grid = std::make_shared<Grid1D>(0, 1, n_cells);
    FemLinearSegment builder(grid);

    // =========== EXPLICIT worker
    {
        // FCT (ncells = 590?)
        ExplicitFct worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.0648930434).margin(1e-5));
    }

    // ============ Crank-Nicolson
    {
        // FCT
        CrankNicolsonFct worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.087307026).margin(1e-5));
    }
}

namespace {

////////////////////////////////////////////////
// FEM-TVD limiter
////////////////////////////////////////////////
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

class TvdAntidiffusion {
public:
    virtual ~TvdAntidiffusion() = default;

    CsrMatrix build_antidiffusion(const std::vector<double>& u) const {
        CsrMatrix F(worker_.fem().stencil());
        const CsrMatrix& L = worker_.low_order_transport();
        const CsrMatrix& D = worker_.transport_diffusion();
        std::vector<double> edges_rij = build_slope_ratio(u);
        for (size_t iedge = 0; iedge < edges_.size(); ++iedge) {
            const DirectedEdge& edge = edges_[iedge];

            double dij = D.vals()[edge.ij_addr];
            double lji = L.vals()[edge.ji_addr];
            double rij = edges_rij[iedge];
            // std::cout << edge.i << ": " << rij << std::endl;

            double fij = std::min(limiter_(rij) * dij, lji);

            F.vals()[edge.ij_addr] += fij;
            F.vals()[edge.ji_addr] += fij;
            F.vals()[edge.ii_addr] -= fij;
            F.vals()[edge.jj_addr] -= fij;
        }
        return F;
    }

protected:
    TvdAntidiffusion(const ATestTransportWorker& worker, std::function<double(double)> limiter)
        : worker_(worker), limiter_(limiter) {
        const CsrMatrix& K = worker_.high_order_transport();
        edges_.clear();
        // build directed edges i->j
        for (auto [i, j, kij]: matrix_iter::ijv(K)) {
            if (i <= j) {
                continue;
            }
            if (kij > 0) {
                std::swap(i, j);
            }
            DirectedEdge edge;
            edge.i = i;
            edge.j = j;
            edge.ij_addr = K.get_address(i, j);
            edge.ji_addr = K.get_address(j, i);
            edge.ii_addr = K.get_address(i, i);
            edge.jj_addr = K.get_address(j, j);
            edges_.push_back(edge);
        }
    }

    const ATestTransportWorker& worker_;
    std::function<double(double)> limiter_;
    struct DirectedEdge {
        size_t i;
        size_t j;
        size_t ij_addr;
        size_t ji_addr;
        size_t ii_addr;
        size_t jj_addr;
    };
    std::vector<DirectedEdge> edges_;

    virtual std::vector<double> build_slope_ratio(const std::vector<double>& u) const = 0;
};

class AlgebraicTvd : public TvdAntidiffusion {
public:
    AlgebraicTvd(const ATestTransportWorker& worker, std::function<double(double)> limiter)
        : TvdAntidiffusion(worker, limiter) {}

    std::vector<double> build_slope_ratio(const std::vector<double>& u) const override {
        const size_t n = u.size();
        const CsrMatrix& K = this->worker_.high_order_transport();
        // Q, P
        std::vector<double> q_plus(n, 0);
        std::vector<double> q_minus(n, 0);
        std::vector<double> p_plus(n, 0);
        std::vector<double> p_minus(n, 0);
        for (auto [i, j, kij]: matrix_iter::ijv(K)) {
            if (i == j) {
                continue;
            }
            double du = u[j] - u[i];
            if (kij > 0) {
                if (du > 0) {
                    q_plus[i] += kij * du;
                } else {
                    q_minus[i] += kij * du;
                }
            } else {
                if (du > 0) {
                    p_minus[i] += kij * du;
                } else {
                    p_plus[i] += kij * du;
                }
            }
        }
        for (size_t i = 0; i < n; ++i) {
            if (std::abs(p_plus[i]) < 1e-12) {
                p_plus[i] = 1e-12;
            }
            if (std::abs(p_minus[i]) < 1e-12) {
                p_minus[i] = 1e-12;
            }
        }
        // Rij
        std::vector<double> ret;
        for (const auto& edge: this->edges_) {
            if (u[edge.i] >= u[edge.j]) {
                ret.push_back(q_plus[edge.i] / p_plus[edge.i]);
            } else {
                ret.push_back(q_minus[edge.i] / p_minus[edge.i]);
            }
        }
        return ret;
    }
};

class InterpolationTvd : public TvdAntidiffusion {
public:
    std::vector<double> build_slope_ratio(const std::vector<double>& u) const override {
        std::vector<double> ret;
        for (size_t iedge = 0; iedge < this->edges_.size(); ++iedge) {
            const auto& edge = this->edges_[iedge];
            double u_ref = u[edge.i];
            double numer = interp_[iedge].apply(u_ref, u);
            double denom = u[edge.j] - u[edge.i];
            if (std::abs(denom) < 1e-12) {
                denom = 1e-12;
            }
            ret.push_back(numer / denom);
        }
        return ret;
    }

protected:
    struct InterpolationData {
        std::vector<std::pair<size_t, double>> weights;
        double apply(double u_ref, const std::vector<double>& u) const {
            double ret = 0;
            for (auto [i, w]: weights) {
                ret += w * (u_ref - u[i]);
            }
            return ret;
        }
    };

    InterpolationTvd(const ATestTransportWorker& worker, std::function<double(double)> limiter)
        : TvdAntidiffusion(worker, limiter) {}

    void set_interpolation_data(std::vector<InterpolationData>&& interp) {
        interp_ = std::move(interp);
    }

private:
    std::vector<InterpolationData> interp_;
};

class ElementBTvd : public InterpolationTvd {
public:
    ElementBTvd(const ATestTransportWorker& worker, std::function<double(double)> limiter)
        : InterpolationTvd(worker, limiter) {
        std::vector<InterpolationData> interp;
        CellFinder cell_finder(worker.grid());
        for (const auto& edge: this->edges_) {
            InterpolationData idata;
            Point pi = worker.grid().point(edge.i);
            Point pj = worker.grid().point(edge.j);
            Point p_upwind = 2 * pi - pj;
            size_t icell = cell_finder(p_upwind);
            if (icell != INVALID_INDEX) {
                Point xi = worker.fem().element(icell).geometry->to_parametric(p_upwind);
                std::vector<size_t> indices = worker.fem().tab_elem_basis(icell);
                std::vector<double> weights = worker.fem().element(icell).basis->value(xi);
                for (size_t i = 0; i < indices.size(); ++i) {
                    idata.weights.push_back(std::make_pair(indices[i], weights[i]));
                }
            }
            interp.push_back(std::move(idata));
        }
        set_interpolation_data(std::move(interp));
    };
};

class GradientTvd : public InterpolationTvd {
public:
    GradientTvd(const ATestTransportWorker& worker, std::function<double(double)> limiter)
        : InterpolationTvd(worker, limiter) {
        const FemAssembler& fem = worker.fem();
        // Compute C matrices
        CsrMatrix cx(fem.stencil());
        CsrMatrix cy(fem.stencil());
        for (size_t ielem = 0; ielem < fem.n_elements(); ++ielem) {
            const FemElement& el = fem.element(ielem);
            std::vector<Vector> c_matrix = element_c_matrix(el);
            auto local_cx = c_matrix | std::views::transform([](Vector v) { return v.x; });
            auto local_cy = c_matrix | std::views::transform([](Vector v) { return v.y; });
            fem.add_to_global_matrix(ielem, std::vector(local_cx.begin(), local_cx.end()), cx.vals());
            fem.add_to_global_matrix(ielem, std::vector(local_cy.begin(), local_cy.end()), cy.vals());
        }

        // build interp
        std::vector<InterpolationData> interp;
        for (const auto& edge: this->edges_) {
            InterpolationData idata;
            Point pi_pj = worker.grid().point(edge.i) - worker.grid().point(edge.j);
            for (auto [k, cx_ik, cy_ik]: matrix_iter::jvv(edge.i, cx, cy)) {
                double w = (pi_pj.x * cx_ik + pi_pj.y * cy_ik) / worker.load_vector()[edge.i];
                idata.weights.push_back(std::make_pair(k, w));
            }
            interp.push_back(std::move(idata));
        }
        set_interpolation_data(std::move(interp));
    }

    std::vector<Vector> element_c_matrix(const FemElement& el) const {
        auto fun = [&el](Point p) -> std::vector<Vector> {
            const size_t n = el.basis->size();
            std::vector<Vector> ret(n * n, Vector{});
            auto val = FemElementValue(&el, p);

            for (size_t ibas = 0; ibas < n; ++ibas) {
                for (size_t jbas = 0; jbas < n; ++jbas) {
                    const size_t k = ibas * n + jbas;
                    ret[k] = val.grad_phi(jbas) * val.phi(ibas) * val.modj();
                }
            }

            return ret;
        };

        return el.quadrature->integrate(fun);
    }
};

template<typename Antidiffusion, Limiter Lim>
class ExplicitTvd : public ATestTransportWorker {
public:
    ExplicitTvd(std::shared_ptr<IGrid> grid, double tau, const IFemBuilder& builder, const ISolution& solution)
        : ATestTransportWorker(grid, tau, builder, solution), antidiffusion_(*this, limiter<Lim>) {}

private:
    Antidiffusion antidiffusion_;

    void impl_step() override {
        CsrMatrix F = antidiffusion_.build_antidiffusion(u_);
        const CsrMatrix& L = low_order_transport();

        std::vector<double> u_new(u_.size(), 0);
        for (size_t i = 0; i < fem_.n_bases(); ++i) {
            double h = load_vector_[i];
            double ku = L.mult_vec(i, u_) - F.mult_vec(i, u_);
            u_new[i] = u_[i] + tau_ / h * ku;
        }
        apply_dirichlet_bc(u_new);
        std::swap(u_, u_new);
    }
};

} // namespace

TEST_CASE("Transport 1D solver, fem-tvd scheme", "[transport1-fem-tvd]") {
    std::cout << std::endl << "--- cfd_test [transport1-fem-tvd] --- " << std::endl;
    // parameters
    const double tend = 0.5;
    const double save_tau = 0.05;
    const double tau = 1e-2;
    const size_t n_cells = 20;

    auto compute = [&](ATestTransportWorker& worker, std::string out_file) {
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

    Solution1D solution;
    auto grid = std::make_shared<Grid1D>(0, 1, n_cells);
    FemLinearSegment builder(grid);

    {
        // UPWIND
        ExplicitTvd<AlgebraicTvd, Limiter::UPWIND> worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.193719381).margin(1e-5));
    }
    {
        // ALGEBRAIC
        ExplicitTvd<AlgebraicTvd, Limiter::MINMOD> worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.1172775615).margin(1e-5));
    }
    {
        // ElementB
        ExplicitTvd<ElementBTvd, Limiter::MINMOD> worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.1172775615).margin(1e-5));
    }
    {
        // Gradient
        ExplicitTvd<GradientTvd, Limiter::MINMOD> worker(grid, tau, builder, solution);
        compute(worker, "transport_fem");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.0857348728).margin(1e-5));
    }
}
