#ifndef __CFD_NUMERIC_INTEGRATION_GAUSSIAN_SQUARE_HPP__
#define __CFD_NUMERIC_INTEGRATION_GAUSSIAN_SQUARE_HPP__

#include "cfd/numeric_integration/quadrature.hpp"
#include <memory>

namespace cfd {

std::shared_ptr<const Quadrature> quadrature_square_gauss1();
std::shared_ptr<const Quadrature> quadrature_square_gauss2();
std::shared_ptr<const Quadrature> quadrature_square_gauss3();
std::shared_ptr<const Quadrature> quadrature_square_gauss4();
std::shared_ptr<const Quadrature> quadrature_square_gauss5();
std::shared_ptr<const Quadrature> quadrature_square_gauss6();

template<int P> std::shared_ptr<const Quadrature> quadrature_square_gauss();

} // namespace cfd
#endif
