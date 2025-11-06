#include "segment_quadrature.hpp"

using namespace cfd;

std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss1() {
    static auto quad = std::make_shared<Quadrature>(std::vector<Point>{{0}}, std::vector<double>{2.0});
    return quad;
}

std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss2() {
    static auto quad = std::make_shared<Quadrature>(std::vector<Point>{{-0.5773502691896257}, {0.5773502691896257}},
                                                    std::vector<double>{1.0, 1.0});

    return quad;
}

std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss3() {
    static auto quad =
        std::make_shared<Quadrature>(std::vector<Point>{{-0.7745966692414834}, {0}, {0.7745966692414834}},
                                     std::vector<double>{0.5555555555555556, 0.8888888888888888, 0.5555555555555556});

    return quad;
}

std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss4() {
    static auto quad = std::make_shared<Quadrature>(
        std::vector<Point>{{-0.8611363115940526}, {-0.3399810435848563}, {0.3399810435848563}, {0.8611363115940526}},
        std::vector<double>{0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538});
    return quad;
}

std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss5() {
    static auto quad =
        std::make_shared<Quadrature>(std::vector<Point>{{0.0000000000000000},
                                                        {-0.5384693101056831},
                                                        {0.5384693101056831},
                                                        {-0.9061798459386640},
                                                        {0.9061798459386640}},
                                     std::vector<double>{0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
                                                         0.2369268850561891, 0.2369268850561891});
    return quad;
}

std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss6() {
    static auto quad =
        std::make_shared<Quadrature>(std::vector<Point>{{0.6612093864662645},
                                                        {-0.6612093864662645},
                                                        {-0.2386191860831969},
                                                        {0.2386191860831969},
                                                        {-0.9324695142031521},
                                                        {0.9324695142031521}},
                                     std::vector<double>{0.3607615730481386, 0.3607615730481386, 0.4679139345726910,
                                                         0.4679139345726910, 0.1713244923791704, 0.1713244923791704});
    return quad;
}

template<int P> std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss() {
    if constexpr (P == 1) {
        return quadrature_segment_gauss1();
    } else if constexpr (P == 2) {
        return quadrature_segment_gauss2();
    } else if constexpr (P == 3) {
        return quadrature_segment_gauss3();
    } else if constexpr (P == 4) {
        return quadrature_segment_gauss4();
    } else if constexpr (P == 5) {
        return quadrature_segment_gauss5();
    } else if constexpr (P == 6) {
        return quadrature_segment_gauss6();
    } else {
        static_assert(false);
    }
}

template std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss<1>();
template std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss<2>();
template std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss<3>();
template std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss<4>();
template std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss<5>();
template std::shared_ptr<const Quadrature> cfd::quadrature_segment_gauss<6>();
