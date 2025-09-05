#ifndef __CFD_NUMERIC_INTEGRATION_GAUSSIAN_SQUARE_HPP__
#define __CFD_NUMERIC_INTEGRATION_GAUSSIAN_SQUARE_HPP__

#include "cfd/numeric_integration/quadrature.hpp"

namespace cfd{

const Quadrature* quadrature_square_gauss1();
const Quadrature* quadrature_square_gauss2();
const Quadrature* quadrature_square_gauss3();
const Quadrature* quadrature_square_gauss4();

}
#endif

