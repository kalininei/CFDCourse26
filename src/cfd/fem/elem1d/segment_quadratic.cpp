#include "segment_quadratic.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Basis
///////////////////////////////////////////////////////////////////////////////
size_t SegmentQuadraticBasis::size() const {
    return 3;
}

std::vector<Point> SegmentQuadraticBasis::parametric_reference_points() const {
    return {Point(-1), Point(1), Point(0)};
}

std::vector<double> SegmentQuadraticBasis::value(Point xi) const {
    double x = xi.x();
    return {(x * x - x) / 2, (x * x + x) / 2, 1 - x * x};
}

std::vector<Vector> SegmentQuadraticBasis::grad(Point xi) const {
    double x = xi.x();
    return {Vector{(2 * x - 1) / 2}, Vector{(2 * x + 1) / 2}, Vector{-(2 * x)}};
}

std::vector<std::array<double, 6>> SegmentQuadraticBasis::upper_hessian(Point xi) const {
    return {{1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0}, {-2, 0, 0, 0, 0, 0}};
}
