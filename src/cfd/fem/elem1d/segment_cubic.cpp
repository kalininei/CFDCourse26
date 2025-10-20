#include "segment_cubic.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Cubic Basis
///////////////////////////////////////////////////////////////////////////////
size_t SegmentCubicBasis::size() const {
    return 4;
}

std::vector<Point> SegmentCubicBasis::parametric_reference_points() const {
    return {Point(-1), Point(1), Point(-1.0 / 3.0), Point(1.0 / 3.0)};
}

std::vector<double> SegmentCubicBasis::value(Point xi) const {
    double x = xi.x();
    return {1.0 / 16.0 * (-(9 * x * x * x) + 9 * x * x + x - 1), 1.0 / 16.0 * (9 * x * x * x + 9 * x * x - x - 1),
            1.0 / 16.0 * (27 * x * x * x - 9 * x * x - 27 * x + 9),
            1.0 / 16.0 * (-(27 * x * x * x) - 9 * x * x + 27 * x + 9)};
}

std::vector<Vector> SegmentCubicBasis::grad(Point xi) const {
    double x = xi.x();
    return {1.0 / 16.0 * Point{-(27 * x * x) + 18 * x + 1}, 1.0 / 16.0 * Point{27 * x * x + 18 * x - 1},
            1.0 / 16.0 * Point{81 * x * x - 18 * x - 27}, 1.0 / 16.0 * Point{-(81 * x * x) - 18 * x + 27}};
}

std::vector<std::array<double, 6>> SegmentCubicBasis::upper_hessian(Point xi) const {
    double x = xi.x();
    return {{1.0 / 16.0 * (-54.0 * x + 18.0), 0, 0, 0, 0, 0},
            {1.0 / 16.0 * (54.0 * x + 18.0), 0, 0, 0, 0, 0},
            {1.0 / 16.0 * (162.0 * x - 18.0), 0, 0, 0, 0, 0},
            {1.0 / 16.0 * (-162.0 * x - 18.0), 0, 0, 0, 0, 0}};
}

///////////////////////////////////////////////////////////////////////////////
// Hermite Basis
///////////////////////////////////////////////////////////////////////////////
SegmentHermiteBasis::SegmentHermiteBasis(std::shared_ptr<const IElementGeometry> geom) : _geom(geom) {}

size_t SegmentHermiteBasis::size() const {
    return 4;
}

std::vector<Point> SegmentHermiteBasis::parametric_reference_points() const {
    return {Point(-1), Point(1), Point(-1), Point(1)};
}

std::vector<double> SegmentHermiteBasis::value(Point xi) const {
    double x = xi.x();
    double C = 1 / _geom->jacobi(xi).modj;
    return {0.25 * (x * x * x - 3 * x + 2), 0.25 * (-x * x * x + 3 * x + 2), 0.25 * ((x * x * x - x * x - x + 1) / C),
            0.25 * ((x * x * x + x * x - x - 1) / C)};
}

std::vector<Vector> SegmentHermiteBasis::grad(Point xi) const {
    double x = xi.x();
    double C = 1 / _geom->jacobi(xi).modj;
    return {0.25 * Point(3 * x * x - 3), 0.25 * Point(3 - 3 * x * x), 0.25 * Point((3 * x * x - 2 * x - 1) / C),
            0.25 * Point((3 * x * x + 2 * x - 1) / C)};
}
