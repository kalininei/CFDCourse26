#include "quadrangle_quadratic.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Basis (9 nodes)
///////////////////////////////////////////////////////////////////////////////

size_t QuadrangleQuadraticBasis::size() const {
    return 9;
}

std::vector<Point> QuadrangleQuadraticBasis::parametric_reference_points() const {
    return {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, -1}, {1, 0}, {0, 1}, {-1, 0}, {0, 0}};
}

std::vector<double> QuadrangleQuadraticBasis::value(Point xi) const {
    auto p0 = [](double x) { return (x * x - x) / 2; };
    auto p1 = [](double x) { return (x * x + x) / 2; };
    auto p2 = [](double x) { return (1 - x * x); };

    double x = xi.x;
    double y = xi.y;
    return {
        p0(x) * p0(y), p1(x) * p0(y), p1(x) * p1(y), p0(x) * p1(y), p2(x) * p0(y),
        p1(x) * p2(y), p2(x) * p1(y), p0(x) * p2(y), p2(x) * p2(y),
    };
}

std::vector<Vector> QuadrangleQuadraticBasis::grad(Point xi) const {
    auto p0 = [](double x) { return (x * x - x) / 2; };
    auto p1 = [](double x) { return (x * x + x) / 2; };
    auto p2 = [](double x) { return (1 - x * x); };
    auto d0 = [](double x) { return (2 * x - 1) / 2; };
    auto d1 = [](double x) { return (2 * x + 1) / 2; };
    auto d2 = [](double x) { return -2 * x; };

    double x = xi.x;
    double y = xi.y;

    return {
        {d0(x) * p0(y), p0(x) * d0(y)}, {d1(x) * p0(y), p1(x) * d0(y)}, {d1(x) * p1(y), p1(x) * d1(y)},
        {d0(x) * p1(y), p0(x) * d1(y)}, {d2(x) * p0(y), p2(x) * d0(y)}, {d1(x) * p2(y), p1(x) * d2(y)},
        {d2(x) * p1(y), p2(x) * d1(y)}, {d0(x) * p2(y), p0(x) * d2(y)}, {d2(x) * p2(y), p2(x) * d2(y)},
    };
}

///////////////////////////////////////////////////////////////////////////////
// Basis (8 nodes)
///////////////////////////////////////////////////////////////////////////////

size_t QuadrangleQuadratic8Basis::size() const {
    return 8;
}

std::vector<Point> QuadrangleQuadratic8Basis::parametric_reference_points() const {
    return {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, -1}, {1, 0}, {0, 1}, {-1, 0}};
}

std::vector<double> QuadrangleQuadratic8Basis::value(Point xi) const {
    double x = xi.x;
    double y = xi.y;
    return {
        0.25 * (-1 + x * y + x * x + y * y - x * x * y - x * y * y),
        0.25 * (-1 - x * y + x * x + y * y - x * x * y + x * y * y),
        0.25 * (-1 + x * y + x * x + y * y + x * x * y + x * y * y),
        0.25 * (-1 - x * y + x * x + y * y + x * x * y - x * y * y),
        0.25 * (2 - 2 * y - 2 * x * x + 2 * x * x * y),
        0.25 * (2 + 2 * x - 2 * y * y - 2 * x * y * y),
        0.25 * (2 + 2 * y - 2 * x * x - 2 * x * x * y),
        0.25 * (2 - 2 * x - 2 * y * y + 2 * x * y * y),
    };
}

std::vector<Vector> QuadrangleQuadratic8Basis::grad(Point xi) const {
    double x = xi.x;
    double y = xi.y;

    return {
        0.25 * Point{y + 2 * x - 2 * x * y - y * y, x + 2 * y - x * x - 2 * x * y},
        0.25 * Point{-y + 2 * x - 2 * x * y + y * y, -x + 2 * y - x * x + 2 * x * y},
        0.25 * Point{y + 2 * x + 2 * x * y + y * y, x + 2 * y + x * x + 2 * x * y},
        0.25 * Point{-y + 2 * x + 2 * x * y - y * y, -x + 2 * y + x * x - 2 * x * y},
        0.25 * Point{-4 * x + 4 * x * y, -2 + 2 * x * x},
        0.25 * Point{2 - 2 * y * y, -4 * y - 4 * x * y},
        0.25 * Point{-4 * x - 4 * x * y, 2 - 2 * x * x},
        0.25 * Point{-2 + 2 * y * y, -4 * y + 4 * x * y},
    };
}
