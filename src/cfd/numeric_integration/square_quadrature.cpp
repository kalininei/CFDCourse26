#include "square_quadrature.hpp"
#include "cfd/numeric_integration/segment_quadrature.hpp"
#include <iostream>

using namespace cfd;

namespace {

Point point_from_segment(const Quadrature *segment_quad, size_t i, size_t j) {
    return Point(segment_quad->points()[i].x(), segment_quad->points()[j].x());
}

double weight_from_segment(const Quadrature *segment_quad, size_t i, size_t j) {
    return segment_quad->weights()[i] * segment_quad->weights()[j];
}

} // namespace

const Quadrature *cfd::quadrature_square_gauss1() {
    static Quadrature quad(
        {
            point_from_segment(quadrature_segment_gauss1(), 0, 0),
        },
        {
            weight_from_segment(quadrature_segment_gauss1(), 0, 0),
        });
    return &quad;
};

const Quadrature *cfd::quadrature_square_gauss2() {
    static Quadrature quad(
        {
            point_from_segment(quadrature_segment_gauss2(), 0, 0),
            point_from_segment(quadrature_segment_gauss2(), 0, 1),
            point_from_segment(quadrature_segment_gauss2(), 1, 0),
            point_from_segment(quadrature_segment_gauss2(), 1, 1),
        },
        {
            weight_from_segment(quadrature_segment_gauss2(), 0, 0),
            weight_from_segment(quadrature_segment_gauss2(), 0, 1),
            weight_from_segment(quadrature_segment_gauss2(), 1, 0),
            weight_from_segment(quadrature_segment_gauss2(), 1, 1),
        });
    return &quad;
};

const Quadrature *cfd::quadrature_square_gauss3() {
    static Quadrature quad(
        {
            point_from_segment(quadrature_segment_gauss3(), 0, 0),
            point_from_segment(quadrature_segment_gauss3(), 0, 1),
            point_from_segment(quadrature_segment_gauss3(), 0, 2),
            point_from_segment(quadrature_segment_gauss3(), 1, 0),
            point_from_segment(quadrature_segment_gauss3(), 1, 1),
            point_from_segment(quadrature_segment_gauss3(), 1, 2),
            point_from_segment(quadrature_segment_gauss3(), 2, 0),
            point_from_segment(quadrature_segment_gauss3(), 2, 1),
            point_from_segment(quadrature_segment_gauss3(), 2, 2),
        },
        {
            weight_from_segment(quadrature_segment_gauss3(), 0, 0),
            weight_from_segment(quadrature_segment_gauss3(), 0, 1),
            weight_from_segment(quadrature_segment_gauss3(), 0, 2),
            weight_from_segment(quadrature_segment_gauss3(), 1, 0),
            weight_from_segment(quadrature_segment_gauss3(), 1, 1),
            weight_from_segment(quadrature_segment_gauss3(), 1, 2),
            weight_from_segment(quadrature_segment_gauss3(), 2, 0),
            weight_from_segment(quadrature_segment_gauss3(), 2, 1),
            weight_from_segment(quadrature_segment_gauss3(), 2, 2),
        });
    return &quad;
};

const Quadrature *cfd::quadrature_square_gauss4() {
    static Quadrature quad(
        {
            point_from_segment(quadrature_segment_gauss4(), 0, 0),
            point_from_segment(quadrature_segment_gauss4(), 0, 1),
            point_from_segment(quadrature_segment_gauss4(), 0, 2),
            point_from_segment(quadrature_segment_gauss4(), 0, 3),
            point_from_segment(quadrature_segment_gauss4(), 1, 0),
            point_from_segment(quadrature_segment_gauss4(), 1, 1),
            point_from_segment(quadrature_segment_gauss4(), 1, 2),
            point_from_segment(quadrature_segment_gauss4(), 1, 3),
            point_from_segment(quadrature_segment_gauss4(), 2, 0),
            point_from_segment(quadrature_segment_gauss4(), 2, 1),
            point_from_segment(quadrature_segment_gauss4(), 2, 2),
            point_from_segment(quadrature_segment_gauss4(), 2, 3),
            point_from_segment(quadrature_segment_gauss4(), 3, 0),
            point_from_segment(quadrature_segment_gauss4(), 3, 1),
            point_from_segment(quadrature_segment_gauss4(), 3, 2),
            point_from_segment(quadrature_segment_gauss4(), 3, 3),
        },
        {
            weight_from_segment(quadrature_segment_gauss4(), 0, 0),
            weight_from_segment(quadrature_segment_gauss4(), 0, 1),
            weight_from_segment(quadrature_segment_gauss4(), 0, 2),
            weight_from_segment(quadrature_segment_gauss4(), 0, 3),
            weight_from_segment(quadrature_segment_gauss4(), 1, 0),
            weight_from_segment(quadrature_segment_gauss4(), 1, 1),
            weight_from_segment(quadrature_segment_gauss4(), 1, 2),
            weight_from_segment(quadrature_segment_gauss4(), 1, 3),
            weight_from_segment(quadrature_segment_gauss4(), 2, 0),
            weight_from_segment(quadrature_segment_gauss4(), 2, 1),
            weight_from_segment(quadrature_segment_gauss4(), 2, 2),
            weight_from_segment(quadrature_segment_gauss4(), 2, 3),
            weight_from_segment(quadrature_segment_gauss4(), 3, 0),
            weight_from_segment(quadrature_segment_gauss4(), 3, 1),
            weight_from_segment(quadrature_segment_gauss4(), 3, 2),
            weight_from_segment(quadrature_segment_gauss4(), 3, 3),
        });
    return &quad;
};
