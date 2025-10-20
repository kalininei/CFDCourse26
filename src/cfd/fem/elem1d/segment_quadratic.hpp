#ifndef __CFD_FEM_SEGMENT_QUADRATIC_HPP__
#define __CFD_FEM_SEGMENT_QUADRATIC_HPP__

#include "cfd/fem/fem_element.hpp"

namespace cfd {

///////////////////////////////////////////////////////////////////////////////
// Basis
///////////////////////////////////////////////////////////////////////////////
class SegmentQuadraticBasis : public IElementBasis {
public:
    size_t size() const override;
    std::vector<Point> parametric_reference_points() const override;
    std::vector<double> value(Point xi) const override;
    std::vector<Vector> grad(Point xi) const override;
    std::vector<std::array<double, 6>> upper_hessian(Point xi) const override;
};

} // namespace cfd
#endif
