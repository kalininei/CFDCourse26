#include "quadrangle_linear.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Geometry
///////////////////////////////////////////////////////////////////////////////
QuadrangleLinearGeometry::QuadrangleLinearGeometry(Point p0, Point p1, Point p2, Point p3)
    : _p0(p0), _p1(p1), _p2(p2), _p3(p3) {
    _c_1 = (p0 + p1 + p2 + p3) / 4;
    _c_xi = (-p0 + p1 + p2 - p3) / 4;
    _c_eta = (-p0 - p1 + p2 + p3) / 4;
    _c_xieta = (p0 - p1 + p2 - p3) / 4;
}

JacobiMatrix QuadrangleLinearGeometry::jacobi(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;

    Point dxi = _c_xi + eta * _c_xieta;
    Point deta = _c_eta + xi * _c_xieta;

    JacobiMatrix jac;
    jac.j11 = dxi.x;
    jac.j12 = deta.x;
    jac.j21 = dxi.y;
    jac.j22 = deta.y;
    jac.j31 = dxi.z;
    jac.j32 = deta.z;

    fill_jacobi_modj_2d(jac);
    return jac;
}

Point QuadrangleLinearGeometry::to_physical(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;
    return _c_1 + xi * _c_xi + eta * _c_eta + xi * eta * _c_xieta;
}

Point QuadrangleLinearGeometry::parametric_center() const {
    return {0, 0};
}

///////////////////////////////////////////////////////////////////////////////
// Basis
///////////////////////////////////////////////////////////////////////////////
size_t QuadrangleLinearBasis::size() const {
    return 4;
}

std::vector<Point> QuadrangleLinearBasis::parametric_reference_points() const {
    return {Point(-1, -1), Point(1, -1), Point(1, 1), Point(-1, 1)};
}

std::vector<double> QuadrangleLinearBasis::value(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;
    return {0.25 * (1 - xi) * (1 - eta), 0.25 * (1 + xi) * (1 - eta), 0.25 * (1 + xi) * (1 + eta),
            0.25 * (1 - xi) * (1 + eta)};
}

std::vector<Vector> QuadrangleLinearBasis::grad(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;
    return {0.25 * Vector(-1 + eta, -1 + xi), 0.25 * Vector(1 - eta, -1 - xi), 0.25 * Vector(1 + eta, 1 + xi),
            0.25 * Vector(-1 - eta, 1 - xi)};
}
