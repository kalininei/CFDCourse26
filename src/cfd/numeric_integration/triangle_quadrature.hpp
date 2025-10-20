#ifndef __CFD_NUMERIC_INTEGRATION_TRIANGLE_QUADRATURE_HPP__
#define __CFD_NUMERIC_INTEGRATION_TRIANGLE_QUADRATURE_HPP__

#include "cfd/numeric_integration/quadrature.hpp"
#include <memory>

namespace cfd {

std::shared_ptr<const Quadrature> quadrature_triangle_gauss1();
std::shared_ptr<const Quadrature> quadrature_triangle_gauss2();
std::shared_ptr<const Quadrature> quadrature_triangle_gauss3();
std::shared_ptr<const Quadrature> quadrature_triangle_gauss4();
std::shared_ptr<const Quadrature> quadrature_triangle_gauss5();
std::shared_ptr<const Quadrature> quadrature_triangle_gauss6();

} // namespace cfd
#endif
