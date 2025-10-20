#ifndef __CFD_FEM_I_ELEMENT_HPP__
#define __CFD_FEM_I_ELEMENT_HPP__

#include "cfd/geom/jacobi.hpp"
#include "cfd/grid/i_grid.hpp"
#include "cfd/mat/densemat.hpp"
#include "cfd/numeric_integration/quadrature.hpp"
#include <functional>

namespace cfd {

///////////////////////////////////////////////////////////////////////////////
// Element Geometry
///////////////////////////////////////////////////////////////////////////////
class IElementGeometry {
public:
    virtual ~IElementGeometry() = default;

    virtual JacobiMatrix jacobi(Point) const = 0;
    virtual Point to_physical(Point) const {
        _THROW_NOT_IMP_;
    }
    virtual Point to_parametric(Point) const {
        _THROW_NOT_IMP_;
    }
    virtual Point parametric_center() const {
        _THROW_NOT_IMP_;
    }
};

///////////////////////////////////////////////////////////////////////////////
// Element Basis
///////////////////////////////////////////////////////////////////////////////
class IElementBasis {
public:
    virtual ~IElementBasis() = default;

    virtual size_t size() const = 0;
    virtual std::vector<Point> parametric_reference_points() const = 0;
    virtual std::vector<double> value(Point) const = 0;
    virtual std::vector<Vector> grad(Point) const = 0;
    virtual std::vector<std::array<double, 6>> upper_hessian(Point) const {
        _THROW_NOT_IMP_;
    }
};

///////////////////////////////////////////////////////////////////////////////
// FemElement
///////////////////////////////////////////////////////////////////////////////
struct FemElement {
    std::shared_ptr<const IElementGeometry> geometry;
    std::shared_ptr<const IElementBasis> basis;
    std::shared_ptr<const Quadrature> quadrature;
};

//////////////////////////////////////////////////////////////////////////////
// FemElementValue
//////////////////////////////////////////////////////////////////////////////
class FemElementValue {
public:
    FemElementValue(const FemElement* elem, Point xi_point);
    ~FemElementValue();

    const IElementGeometry* geom;
    const IElementBasis* basis;
    Point xi;

    double phi(size_t i) const;      // phi_i at (xi)
    Vector grad_phi(size_t i) const; // dphi_i/dx,dy,dz at xi
    double laplace(size_t i) const;  // div(grad(phi_i)) at xi
    double modj() const;
    const JacobiMatrix* jacobi() const;
    double interpolate(const std::vector<double>& f) const;
    Vector interpolate(const std::vector<Vector>& f) const;
    double divergence(const std::vector<Vector>& f) const;

private:
    struct Cache;
    std::unique_ptr<Cache> pcache_;
};

} // namespace cfd
#endif
