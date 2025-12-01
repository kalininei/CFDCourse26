#include "cfd/fem/elem1d/segment_linear.hpp"
#include "cfd/fem/elem2d/triangle_linear.hpp"
#include "cfd/fem/fem_assembler.hpp"
#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd/mat/csrmat.hpp"
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

///////////////////////////////////////////////////////////
// BasicWorker
///////////////////////////////////////////////////////////
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

        // velocity
        for (size_t ibas = 0; ibas < fem_.n_bases(); ++ibas) {
            Point p = fem_.reference_point(ibas);
            velocity_.push_back(solution_.velocity(p));
        }
        // initial solution
        for (size_t ibas = 0; ibas < u_.size(); ++ibas) {
            Point p = fem_.reference_point(ibas);
            u_[ibas] = init_solution(p);
        }
        // load vector
        load_vector_.resize(fem_.n_bases(), 0);
        for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
            std::vector<double> local = element_load_vector(ielem);
            fem_.add_to_global_vector(ielem, local, load_vector_);
        }
        // stabilized mass
        stabilized_mass_matrix_.set_stencil(fem_.stencil());
        for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
            std::vector<double> local = element_stabilized_mass_matrix(ielem);
            fem_.add_to_global_matrix(ielem, local, stabilized_mass_matrix_.vals());
        }
        // stabilized transport
        stabilized_transport_matrix_.set_stencil(fem_.stencil());
        for (size_t ielem = 0; ielem < fem_.n_elements(); ++ielem) {
            std::vector<double> local = element_stabilized_transport_matrix(ielem);
            fem_.add_to_global_matrix(ielem, local, stabilized_transport_matrix_.vals());
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

    const IGrid& grid() const {
        return *grid_;
    }

protected:
    std::shared_ptr<IGrid> grid_;
    FemAssembler fem_;
    std::vector<double> load_vector_;
    CsrMatrix stabilized_mass_matrix_;
    CsrMatrix stabilized_transport_matrix_;

    const double tau_;
    std::vector<double> u_;
    std::vector<Vector> velocity_;
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
    void apply_dirichlet_bc(CsrMatrix& m) const {
        for (size_t ibas: boundary_bases()) {
            m.set_unit_row(ibas);
        }
    }

private:
    virtual void impl_step() = 0;

    double stab(size_t ielem, size_t ibas, const FemElementValue& val) const {
        std::vector<Vector> vel_vec = fem_.local_vector(ielem, velocity_);
        Vector velocity = val.interpolate(vel_vec);
        double abs_velocity = vector_abs(velocity);

        double h = grid_->cell_size(ielem);
        double eps = abs_velocity * h / 2;

        Vector grad_phi = val.grad_phi(ibas);
        return eps / abs_velocity / abs_velocity * dot_product(velocity, grad_phi);
    }

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

    std::vector<double> element_stabilized_mass_matrix(size_t ielem) const {
        const FemElement& el = fem_.element(ielem);

        auto fun = [&](Point p) -> std::vector<double> {
            const size_t n = el.basis->size();
            std::vector<double> ret(n * n, 0.0);
            auto val = FemElementValue(&el, p);

            for (size_t ibas = 0; ibas < n; ++ibas) {
                for (size_t jbas = 0; jbas < n; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    ret[k1] = val.phi(jbas) * (val.phi(ibas) + stab(ielem, ibas, val)) * val.modj();
                }
            }

            return ret;
        };

        return el.quadrature->integrate(fun);
    }

    std::vector<double> element_stabilized_transport_matrix(size_t ielem) const {
        const FemElement& el = fem_.element(ielem);

        auto fun = [&](Point p) -> std::vector<double> {
            const size_t n = el.basis->size();
            const std::vector<Vector> element_velocity = fem_.local_vector(ielem, velocity_);

            std::vector<double> ret(n * n, 0.0);
            auto val = FemElementValue(&el, p);

            for (size_t ibas = 0; ibas < n; ++ibas) {
                for (size_t jbas = 0; jbas < n; ++jbas) {
                    const size_t k1 = ibas * n + jbas;
                    Vector vel = val.interpolate(element_velocity);
                    Vector g2 = val.grad_phi(jbas);
                    ret[k1] = -dot_product(vel, g2) * (val.phi(ibas) + stab(ielem, ibas, val)) * val.modj();
                }
            }
            return ret;
        };
        return el.quadrature->integrate(fun);
    }
};

/////////////////////////////////////////////////////////
// ThetaSupg
/////////////////////////////////////////////////////////
class ThetaSupg : public ATestTransportWorker {
public:
    ThetaSupg(std::shared_ptr<IGrid> grid, double tau, double theta, const IFemBuilder& builder,
              const ISolution& solution)
        : ATestTransportWorker(grid, tau, builder, solution), solver_(1000, 1e-12) {

        lhs_ = fem_.zero_matrix();
        rhs_ = fem_.zero_matrix();
        for (auto [m_ij, k_ij, lhs_ij, rhs_ij]:
             matrix_iter::v(stabilized_mass_matrix_, stabilized_transport_matrix_, lhs_, rhs_)) {
            lhs_ij = m_ij - tau * theta * k_ij;
            rhs_ij = m_ij + tau * (1 - theta) * k_ij;
        }

        apply_dirichlet_bc(lhs_);
        solver_.set_matrix(lhs_);
    }

private:
    CsrMatrix lhs_;
    CsrMatrix rhs_;
    AmgcMatrixSolver solver_;

    void impl_step() override {
        std::vector<double> f = rhs_.mult_vec(u_);
        apply_dirichlet_bc(f);
        solver_.solve(f, u_);
    }
};

} // namespace

TEST_CASE("SUPG 1d", "[transport1-supg]") {
    std::cout << std::endl << "--- cfd_test [transport1-supg] --- " << std::endl;
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
        // Crank-Nicolson
        ThetaSupg worker(grid, tau, 0.5, builder, solution);
        compute(worker, "transport_supg");
        double norm = worker.compute_norm2();
        std::cout << n_cells << " " << norm << std::endl;
        CHECK(norm == Approx(0.0473896818).margin(1e-5));
    }
}
