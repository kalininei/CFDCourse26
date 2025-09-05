#ifndef __CFD_NUMERIC_INTEGRATION_GAUSSIAN_SEGMENT_HPP__
#define __CFD_NUMERIC_INTEGRATION_GAUSSIAN_SEGMENT_HPP__

#include "cfd/numeric_integration/quadrature.hpp"

namespace cfd{

const Quadrature* quadrature_segment_gauss1();
const Quadrature* quadrature_segment_gauss2();
const Quadrature* quadrature_segment_gauss3();
const Quadrature* quadrature_segment_gauss4();
const Quadrature* quadrature_segment_gauss5();
const Quadrature* quadrature_segment_gauss6();

}
#endif
