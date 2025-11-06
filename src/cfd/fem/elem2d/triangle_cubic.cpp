#include "triangle_cubic.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Basis, 10 node cubic
///////////////////////////////////////////////////////////////////////////////
size_t TriangleCubicBasis::size() const {
    return 10;
}

std::vector<Point> TriangleCubicBasis::parametric_reference_points() const {
    return {{0, 0},
            {1, 0},
            {0, 1},
            {1.0 / 3.0, 0},
            {2.0 / 3.0, 0},
            {2.0 / 3.0, 1.0 / 3.0},
            {1.0 / 3.0, 2.0 / 3.0},
            {0, 2.0 / 3.0},
            {0, 1.0 / 3.0},
            {1.0 / 3.0, 1.0 / 3.0}};
}

std::vector<double> TriangleCubicBasis::value(Point xi) const {
    double x = xi.x;
    double y = xi.y;
    return {
        0.5 * (-9 * y * y * y + (18 - 27 * x) * y * y + (-27 * x * x + 36 * x - 11) * y - 9 * x * x * x + 18 * x * x -
               11 * x + 2),
        0.5 * (9 * x * x * x - 9 * x * x + 2 * x),
        0.5 * (9 * y * y * y - 9 * y * y + 2 * y),
        0.5 * (27 * x * y * y + (54 * x * x - 45 * x) * y + 27 * x * x * x - 45 * x * x + 18 * x),
        0.5 * ((9 * x - 27 * x * x) * y - 27 * x * x * x + 36 * x * x - 9 * x),
        0.5 * ((27 * x * x - 9 * x) * y),
        0.5 * (27 * x * y * y - 9 * x * y),
        0.5 * (-27 * y * y * y + (36 - 27 * x) * y * y + (9 * x - 9) * y),
        0.5 * (27 * y * y * y + (54 * x - 45) * y * y + (27 * x * x - 45 * x + 18) * y),
        0.5 * ((54 * x - 54 * x * x) * y - 54 * x * y * y),
    };
}

std::vector<Vector> TriangleCubicBasis::grad(Point xi) const {
    double x = xi.x;
    double y = xi.y;
    return {0.5 * Vector{-27 * y * y + (36 - 54 * x) * y - 27 * x * x + 36 * x - 11,
                         -27 * y * y + (36 - 54 * x) * y - 27 * x * x + 36 * x - 11},
            0.5 * Vector{27 * x * x - 18 * x + 2, 0},
            0.5 * Vector{0, 27 * y * y - 18 * y + 2},
            0.5 * Vector{27 * y * y + (108 * x - 45) * y + 81 * x * x - 90 * x + 18, 54 * x * y + 54 * x * x - 45 * x},
            0.5 * Vector{(9 - 54 * x) * y - 81 * x * x + 72 * x - 9, 9 * x - 27 * x * x},
            0.5 * Vector{(54 * x - 9) * y, 27 * x * x - 9 * x},
            0.5 * Vector{27 * y * y - 9 * y, 54 * x * y - 9 * x},
            0.5 * Vector{9 * y - 27 * y * y, -81 * y * y + (72 - 54 * x) * y + 9 * x - 9},
            0.5 * Vector{54 * y * y + (54 * x - 45) * y, 81 * y * y + (108 * x - 90) * y + 27 * x * x - 45 * x + 18},
            0.5 * Vector{(54 - 108 * x) * y - 54 * y * y, -108 * x * y - 54 * x * x + 54 * x}};
}
