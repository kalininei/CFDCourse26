#include "square_quadrature.hpp"
#include "cfd/numeric_integration/segment_quadrature.hpp"
#include <iostream>

using namespace cfd;

namespace {

Point point_from_segment(std::shared_ptr<const Quadrature> segment_quad, size_t i, size_t j) {
    return Point(segment_quad->points()[i].x, segment_quad->points()[j].x);
}

double weight_from_segment(std::shared_ptr<const Quadrature> segment_quad, size_t i, size_t j) {
    return segment_quad->weights()[i] * segment_quad->weights()[j];
}

} // namespace

std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss1() {
    static auto quad = std::make_shared<Quadrature>(
        std::vector<Point>{
            point_from_segment(quadrature_segment_gauss1(), 0, 0),
        },
        std::vector<double>{
            weight_from_segment(quadrature_segment_gauss1(), 0, 0),
        });
    return quad;
};

std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss2() {
    static auto quad = std::make_shared<Quadrature>(
        std::vector<Point>{
            point_from_segment(quadrature_segment_gauss2(), 0, 0),
            point_from_segment(quadrature_segment_gauss2(), 0, 1),
            point_from_segment(quadrature_segment_gauss2(), 1, 0),
            point_from_segment(quadrature_segment_gauss2(), 1, 1),
        },
        std::vector<double>{
            weight_from_segment(quadrature_segment_gauss2(), 0, 0),
            weight_from_segment(quadrature_segment_gauss2(), 0, 1),
            weight_from_segment(quadrature_segment_gauss2(), 1, 0),
            weight_from_segment(quadrature_segment_gauss2(), 1, 1),
        });
    return quad;
};

std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss3() {
    static auto quad = std::make_shared<Quadrature>(
        std::vector<Point>{
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
        std::vector<double>{
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
    return quad;
};

std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss4() {
    static auto quad = std::make_shared<Quadrature>(
        std::vector<Point>{
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
        std::vector<double>{
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
    return quad;
};

std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss5() {
    static auto quad = std::make_shared<Quadrature>(
        std::vector<Point>{
            point_from_segment(quadrature_segment_gauss5(), 0, 0),
            point_from_segment(quadrature_segment_gauss5(), 0, 1),
            point_from_segment(quadrature_segment_gauss5(), 0, 2),
            point_from_segment(quadrature_segment_gauss5(), 0, 3),
            point_from_segment(quadrature_segment_gauss5(), 0, 4),
            point_from_segment(quadrature_segment_gauss5(), 1, 0),
            point_from_segment(quadrature_segment_gauss5(), 1, 1),
            point_from_segment(quadrature_segment_gauss5(), 1, 2),
            point_from_segment(quadrature_segment_gauss5(), 1, 3),
            point_from_segment(quadrature_segment_gauss5(), 1, 4),
            point_from_segment(quadrature_segment_gauss5(), 2, 0),
            point_from_segment(quadrature_segment_gauss5(), 2, 1),
            point_from_segment(quadrature_segment_gauss5(), 2, 2),
            point_from_segment(quadrature_segment_gauss5(), 2, 3),
            point_from_segment(quadrature_segment_gauss5(), 2, 4),
            point_from_segment(quadrature_segment_gauss5(), 3, 0),
            point_from_segment(quadrature_segment_gauss5(), 3, 1),
            point_from_segment(quadrature_segment_gauss5(), 3, 2),
            point_from_segment(quadrature_segment_gauss5(), 3, 3),
            point_from_segment(quadrature_segment_gauss5(), 3, 4),
            point_from_segment(quadrature_segment_gauss5(), 4, 0),
            point_from_segment(quadrature_segment_gauss5(), 4, 1),
            point_from_segment(quadrature_segment_gauss5(), 4, 2),
            point_from_segment(quadrature_segment_gauss5(), 4, 3),
            point_from_segment(quadrature_segment_gauss5(), 4, 4),
        },
        std::vector<double>{
            weight_from_segment(quadrature_segment_gauss5(), 0, 0),
            weight_from_segment(quadrature_segment_gauss5(), 0, 1),
            weight_from_segment(quadrature_segment_gauss5(), 0, 2),
            weight_from_segment(quadrature_segment_gauss5(), 0, 3),
            weight_from_segment(quadrature_segment_gauss5(), 0, 4),
            weight_from_segment(quadrature_segment_gauss5(), 1, 0),
            weight_from_segment(quadrature_segment_gauss5(), 1, 1),
            weight_from_segment(quadrature_segment_gauss5(), 1, 2),
            weight_from_segment(quadrature_segment_gauss5(), 1, 3),
            weight_from_segment(quadrature_segment_gauss5(), 1, 4),
            weight_from_segment(quadrature_segment_gauss5(), 2, 0),
            weight_from_segment(quadrature_segment_gauss5(), 2, 1),
            weight_from_segment(quadrature_segment_gauss5(), 2, 2),
            weight_from_segment(quadrature_segment_gauss5(), 2, 3),
            weight_from_segment(quadrature_segment_gauss5(), 2, 4),
            weight_from_segment(quadrature_segment_gauss5(), 3, 0),
            weight_from_segment(quadrature_segment_gauss5(), 3, 1),
            weight_from_segment(quadrature_segment_gauss5(), 3, 2),
            weight_from_segment(quadrature_segment_gauss5(), 3, 3),
            weight_from_segment(quadrature_segment_gauss5(), 3, 4),
            weight_from_segment(quadrature_segment_gauss5(), 4, 0),
            weight_from_segment(quadrature_segment_gauss5(), 4, 1),
            weight_from_segment(quadrature_segment_gauss5(), 4, 2),
            weight_from_segment(quadrature_segment_gauss5(), 4, 3),
            weight_from_segment(quadrature_segment_gauss5(), 4, 4),
        });
    return quad;
};

std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss6() {
    static auto quad = std::make_shared<Quadrature>(
        std::vector<Point>{
            point_from_segment(quadrature_segment_gauss6(), 0, 0),
            point_from_segment(quadrature_segment_gauss6(), 0, 1),
            point_from_segment(quadrature_segment_gauss6(), 0, 2),
            point_from_segment(quadrature_segment_gauss6(), 0, 3),
            point_from_segment(quadrature_segment_gauss6(), 0, 4),
            point_from_segment(quadrature_segment_gauss6(), 0, 5),
            point_from_segment(quadrature_segment_gauss6(), 1, 0),
            point_from_segment(quadrature_segment_gauss6(), 1, 1),
            point_from_segment(quadrature_segment_gauss6(), 1, 2),
            point_from_segment(quadrature_segment_gauss6(), 1, 3),
            point_from_segment(quadrature_segment_gauss6(), 1, 4),
            point_from_segment(quadrature_segment_gauss6(), 1, 5),
            point_from_segment(quadrature_segment_gauss6(), 2, 0),
            point_from_segment(quadrature_segment_gauss6(), 2, 1),
            point_from_segment(quadrature_segment_gauss6(), 2, 2),
            point_from_segment(quadrature_segment_gauss6(), 2, 3),
            point_from_segment(quadrature_segment_gauss6(), 2, 4),
            point_from_segment(quadrature_segment_gauss6(), 2, 5),
            point_from_segment(quadrature_segment_gauss6(), 3, 0),
            point_from_segment(quadrature_segment_gauss6(), 3, 1),
            point_from_segment(quadrature_segment_gauss6(), 3, 2),
            point_from_segment(quadrature_segment_gauss6(), 3, 3),
            point_from_segment(quadrature_segment_gauss6(), 3, 4),
            point_from_segment(quadrature_segment_gauss6(), 3, 5),
            point_from_segment(quadrature_segment_gauss6(), 4, 0),
            point_from_segment(quadrature_segment_gauss6(), 4, 1),
            point_from_segment(quadrature_segment_gauss6(), 4, 2),
            point_from_segment(quadrature_segment_gauss6(), 4, 3),
            point_from_segment(quadrature_segment_gauss6(), 4, 4),
            point_from_segment(quadrature_segment_gauss6(), 4, 5),
            point_from_segment(quadrature_segment_gauss6(), 5, 0),
            point_from_segment(quadrature_segment_gauss6(), 5, 1),
            point_from_segment(quadrature_segment_gauss6(), 5, 2),
            point_from_segment(quadrature_segment_gauss6(), 5, 3),
            point_from_segment(quadrature_segment_gauss6(), 5, 4),
            point_from_segment(quadrature_segment_gauss6(), 5, 5),
        },
        std::vector<double>{
            weight_from_segment(quadrature_segment_gauss6(), 0, 0),
            weight_from_segment(quadrature_segment_gauss6(), 0, 1),
            weight_from_segment(quadrature_segment_gauss6(), 0, 2),
            weight_from_segment(quadrature_segment_gauss6(), 0, 3),
            weight_from_segment(quadrature_segment_gauss6(), 0, 4),
            weight_from_segment(quadrature_segment_gauss6(), 0, 5),
            weight_from_segment(quadrature_segment_gauss6(), 1, 0),
            weight_from_segment(quadrature_segment_gauss6(), 1, 1),
            weight_from_segment(quadrature_segment_gauss6(), 1, 2),
            weight_from_segment(quadrature_segment_gauss6(), 1, 3),
            weight_from_segment(quadrature_segment_gauss6(), 1, 4),
            weight_from_segment(quadrature_segment_gauss6(), 1, 5),
            weight_from_segment(quadrature_segment_gauss6(), 2, 0),
            weight_from_segment(quadrature_segment_gauss6(), 2, 1),
            weight_from_segment(quadrature_segment_gauss6(), 2, 2),
            weight_from_segment(quadrature_segment_gauss6(), 2, 3),
            weight_from_segment(quadrature_segment_gauss6(), 2, 4),
            weight_from_segment(quadrature_segment_gauss6(), 2, 5),
            weight_from_segment(quadrature_segment_gauss6(), 3, 0),
            weight_from_segment(quadrature_segment_gauss6(), 3, 1),
            weight_from_segment(quadrature_segment_gauss6(), 3, 2),
            weight_from_segment(quadrature_segment_gauss6(), 3, 3),
            weight_from_segment(quadrature_segment_gauss6(), 3, 4),
            weight_from_segment(quadrature_segment_gauss6(), 3, 5),
            weight_from_segment(quadrature_segment_gauss6(), 4, 0),
            weight_from_segment(quadrature_segment_gauss6(), 4, 1),
            weight_from_segment(quadrature_segment_gauss6(), 4, 2),
            weight_from_segment(quadrature_segment_gauss6(), 4, 3),
            weight_from_segment(quadrature_segment_gauss6(), 4, 4),
            weight_from_segment(quadrature_segment_gauss6(), 4, 5),
            weight_from_segment(quadrature_segment_gauss6(), 5, 0),
            weight_from_segment(quadrature_segment_gauss6(), 5, 1),
            weight_from_segment(quadrature_segment_gauss6(), 5, 2),
            weight_from_segment(quadrature_segment_gauss6(), 5, 3),
            weight_from_segment(quadrature_segment_gauss6(), 5, 4),
            weight_from_segment(quadrature_segment_gauss6(), 5, 5),
        });
    return quad;
};

template<int P>
std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss() {
    if constexpr (P == 1) {
        return quadrature_square_gauss1();
    } else if constexpr (P == 2) {
        return quadrature_square_gauss2();
    } else if constexpr (P == 3) {
        return quadrature_square_gauss3();
    } else if constexpr (P == 4) {
        return quadrature_square_gauss4();
    } else if constexpr (P == 5) {
        return quadrature_square_gauss5();
    } else if constexpr (P == 6) {
        return quadrature_square_gauss6();
    } else {
        static_assert(false);
    }
}

template std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss<1>();
template std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss<2>();
template std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss<3>();
template std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss<4>();
template std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss<5>();
template std::shared_ptr<const Quadrature> cfd::quadrature_square_gauss<6>();
