#include "segment_linear.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Geometry
///////////////////////////////////////////////////////////////////////////////
SegmentLinearGeometry::SegmentLinearGeometry(Point p0, Point p1) : _p0(p0), _p1(p1) {
    double modj = vector_abs(p1 - p0) / 2;
    _jac.j11 = modj;
    _jac.modj = modj;
}

JacobiMatrix SegmentLinearGeometry::jacobi(Point) const {
    return _jac;
}

Point SegmentLinearGeometry::to_physical(Point xi) const {
    double t = (xi.x + 1) / 2.0;
    return (1 - t) * _p0 + t * _p1;
}

Point SegmentLinearGeometry::to_parametric(Point p) const {
    double t = (p.x - _p0.x) / (_p1.x - _p0.x);
    return Point{2 * t - 1};
}

Point SegmentLinearGeometry::parametric_center() const {
    return {0};
}

///////////////////////////////////////////////////////////////////////////////
// Basis
///////////////////////////////////////////////////////////////////////////////
size_t SegmentLinearBasis::size() const {
    return 2;
}

std::vector<Point> SegmentLinearBasis::parametric_reference_points() const {
    return {Point(-1), Point(1)};
}

std::vector<double> SegmentLinearBasis::value(Point xi_) const {
    double xi = xi_.x;
    return {(1 - xi) / 2, (1 + xi) / 2};
}

std::vector<Vector> SegmentLinearBasis::grad(Point) const {
    return {Vector{-0.5}, Vector{0.5}};
}

std::vector<std::array<double, 6>> SegmentLinearBasis::upper_hessian(Point) const {
    return {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
}
