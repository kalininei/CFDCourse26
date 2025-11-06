#include "fem_boundary.hpp"

using namespace cfd;

namespace {

class BoundaryElementGeometry0D : public IElementGeometry {
public:
    BoundaryElementGeometry0D() {}

    JacobiMatrix jacobi(Point) const override {
        const double modj = 1.0;
        return {modj, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, modj};
    }
};

class BoundaryElementGeometry1D : public IElementGeometry {
public:
    BoundaryElementGeometry1D(std::shared_ptr<const IElementGeometry> internal_geom, Vector jtau,
                              std::function<Point(Point)> xi_from_tau)
        : internal_geom_(internal_geom), jtau_(jtau), xi_from_tau_(xi_from_tau) {}

    JacobiMatrix jacobi(Point tau) const override {
        const Point xi = xi_from_tau_(tau);
        const JacobiMatrix internal_j = internal_geom_->jacobi(xi);
        const Vector v = {internal_j.j11 * jtau_.x + internal_j.j12 * jtau_.y,
                          internal_j.j21 * jtau_.x + internal_j.j22 * jtau_.y};
        const double modj = vector_abs(v);
        return {modj, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, modj};
    }

private:
    std::shared_ptr<const IElementGeometry> internal_geom_;
    Vector jtau_;
    std::function<Point(Point)> xi_from_tau_;
};

class BoundaryElementBasis : public IElementBasis {
public:
    BoundaryElementBasis(std::shared_ptr<const IElementBasis> internal_basis, std::function<Point(Point)> xi_from_tau)
        : internal_(internal_basis), xi_from_tau_(xi_from_tau) {}

    size_t size() const override {
        return internal_->size();
    }

    std::vector<Point> parametric_reference_points() const override {
        _THROW_NOT_IMP_;
    }

    std::vector<double> value(Point p) const override {
        return internal_->value(xi_from_tau_(p));
    }
    std::vector<Vector> grad(Point) const override {
        _THROW_NOT_IMP_;
    };

private:
    std::shared_ptr<const IElementBasis> internal_;
    std::function<Point(Point)> xi_from_tau_;
};

} // namespace

FemElement cfd::build_boundary_element(const FemElement& internal_element, std::shared_ptr<const Quadrature> quad,
                                       std::function<Point(Point)> xi_from_tau, Vector Jtau, Vector Jsigma) {

    std::shared_ptr<IElementGeometry> geom;
    std::shared_ptr<IElementBasis> basis;

    if (vector_abs(Jtau) < 1e-16 && vector_abs(Jsigma) < 1e-16) {
        // boundary of 1d grid
        geom = std::make_shared<BoundaryElementGeometry0D>();
    } else if (vector_abs(Jsigma) < 1e-16) {
        // boundary of 2d grid
        geom = std::make_shared<BoundaryElementGeometry1D>(internal_element.geometry, Jtau, xi_from_tau);
    }

    basis = std::make_shared<BoundaryElementBasis>(internal_element.basis, xi_from_tau);

    return FemElement{geom, basis, quad};
}
