#include "triangle_linear.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Geometry
///////////////////////////////////////////////////////////////////////////////
TriangleLinearGeometry::TriangleLinearGeometry(Point p0, Point p1, Point p2) : _p0(p0), _p1(p1), _p2(p2) {
    _jac.modj = vector_abs(cross_product(p1 - p0, p2 - p0));
    _jac.j11 = p1.x - p0.x;
    _jac.j12 = p2.x - p0.x;
    _jac.j21 = p1.y - p0.y;
    _jac.j22 = p2.y - p0.y;
    _jac.j31 = p1.z - p0.z;
    _jac.j32 = p2.z - p0.z;
}

JacobiMatrix TriangleLinearGeometry::jacobi(Point) const {
    return _jac;
}

Point TriangleLinearGeometry::to_physical(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;
    return (1 - xi - eta) * _p0 + xi * _p1 + eta * _p2;
}

Point TriangleLinearGeometry::to_parametric(Point) const {
    _THROW_NOT_IMP_;
}

Point TriangleLinearGeometry::parametric_center() const {
    return {1.0 / 3.0, 1.0 / 3.0};
}

///////////////////////////////////////////////////////////////////////////////
// Basis
///////////////////////////////////////////////////////////////////////////////
size_t TriangleLinearBasis::size() const {
    return 3;
}

std::vector<Point> TriangleLinearBasis::parametric_reference_points() const {
    return {Point(0, 0), Point(1, 0), Point(0, 1)};
}

std::vector<double> TriangleLinearBasis::value(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;
    return {1 - xi - eta, xi, eta};
}

std::vector<Vector> TriangleLinearBasis::grad(Point) const {
    return {Vector(-1, -1), Vector(1, 0), Vector(0, 1)};
}
