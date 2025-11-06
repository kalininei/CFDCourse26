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
// QuadrangleLinearBasis
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

///////////////////////////////////////////////////////////////////////////////
// Quadrangle_Sqr0_Linear123
///////////////////////////////////////////////////////////////////////////////
size_t Quadrangle_Sqr0_Linear123::size() const {
    return 5;
}

std::vector<Point> Quadrangle_Sqr0_Linear123::parametric_reference_points() const {
    return {Point(-1, -1), Point(1, -1), Point(1, 1), Point(-1, 1), Point(0, -1)};
}

std::vector<double> Quadrangle_Sqr0_Linear123::value(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;
    double n4 = 0.5 * (xi - 1) * (xi + 1) * (eta - 1);
    std::vector<double> n_prime = linear_.value(xi_);

    return {n_prime[0] - 0.5 * n4, n_prime[1] - 0.5 * n4, n_prime[2], n_prime[3], n4};
}

std::vector<Vector> Quadrangle_Sqr0_Linear123::grad(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;
    Vector grad4(xi * (eta - 1), (xi - 1) * (xi + 1) / 2.0);
    std::vector<Vector> grad_prime = linear_.grad(xi_);

    return {grad_prime[0] - 0.5 * grad4, grad_prime[1] - 0.5 * grad4, grad_prime[2], grad_prime[3], grad4};
}

///////////////////////////////////////////////////////////////////////////////
// Quadrangle_Cub0_Linear123
///////////////////////////////////////////////////////////////////////////////
size_t Quadrangle_Cub0_Linear123::size() const {
    return 6;
}

std::vector<Point> Quadrangle_Cub0_Linear123::parametric_reference_points() const {
    return {Point(-1, -1), Point(1, -1), Point(1, 1), Point(-1, 1), Point(-1.0 / 3.0, -1), Point(1.0 / 3.0, -1)};
}

std::vector<double> Quadrangle_Cub0_Linear123::value(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;

    double n4 = (1.0 - 3.0 * xi) * (xi - 1.0) * (xi + 1.0) * (eta - 1.0) * 9.0 / 32.0;
    double n5 = (1.0 + 3.0 * xi) * (xi - 1.0) * (xi + 1.0) * (eta - 1.0) * 9.0 / 32.0;

    std::vector<double> n_prime = linear_.value(xi_);

    return {
        n_prime[0] - 2.0 / 3.0 * n4 - 1.0 / 3.0 * n5, //
        n_prime[1] - 1.0 / 3.0 * n4 - 2.0 / 3.0 * n5, //
        n_prime[2],                                   //
        n_prime[3],                                   //
        n4,                                           //
        n5                                            //
    };
}

std::vector<Vector> Quadrangle_Cub0_Linear123::grad(Point xi_) const {
    double xi = xi_.x;
    double eta = xi_.y;

    Vector grad4{
        9.0 / 32.0 * (eta - 1.0) * (-3.0 * (xi * xi - 1.0) + 2 * xi * (1.0 - 3.0 * xi)), //
        (1.0 - 3.0 * xi) * (xi - 1.0) * (xi + 1.0) * 9.0 / 32.0                          //
    };
    Vector grad5{
        9.0 / 32.0 * (eta - 1.0) * (+3.0 * (xi * xi - 1.0) + 2 * xi * (1.0 + 3.0 * xi)), //
        (1.0 + 3.0 * xi) * (xi - 1.0) * (xi + 1.0) * 9.0 / 32.0                          //
    };
    std::vector<Vector> grad_prime = linear_.grad(xi_);

    return {
        grad_prime[0] - 2.0 / 3.0 * grad4 - 1.0 / 3.0 * grad5, //
        grad_prime[1] - 1.0 / 3.0 * grad4 - 2.0 / 3.0 * grad5, //
        grad_prime[2],                                         //
        grad_prime[3],                                         //
        grad4,                                                 //
        grad5                                                  //
    };
}
