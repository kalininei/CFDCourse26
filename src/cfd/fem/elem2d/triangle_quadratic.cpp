#include "triangle_quadratic.hpp"

using namespace cfd;

size_t TriangleQuadraticBasis::size() const {
    return 6;
}

std::vector<Point> TriangleQuadraticBasis::parametric_reference_points() const {
    return {{0, 0}, {1, 0}, {0, 1}, {0.5, 0}, {0.5, 0.5}, {0, 0.5}};
}

std::vector<double> TriangleQuadraticBasis::value(Point xi) const {
    double x = xi.x;
    double y = xi.y;
    return {
        1 - 3 * x - 3 * y + 4 * x * y + 2 * x * x + 2 * y * y,
        -x + 2 * x * x,
        -y + 2 * y * y,
        4 * x - 4 * x * y - 4 * x * x,
        4 * x * y,
        4 * y - 4 * x * y - 4 * y * y,
    };
}

std::vector<Vector> TriangleQuadraticBasis::grad(Point xi) const {
    double x = xi.x;
    double y = xi.y;

    return {{-3 + 4 * y + 4 * x, -3 + 4 * x + 4 * y},
            {-1 + 4 * x, 0},
            {0, -1 + 4 * y},
            {4 - 4 * y - 8 * x, -4 * x},
            {4 * y, 4 * x},
            {-4 * y, 4 - 4 * x - 8 * y}};
}
