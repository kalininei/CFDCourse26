#ifndef __CFD_FEM_BOUNDARY_HPP__
#define __CFD_FEM_BOUNDARY_HPP__

#include "cfd/fem/fem_element.hpp"

namespace cfd {

FemElement build_boundary_element(const FemElement& internal_element, std::shared_ptr<const Quadrature> quadrature,
                                  std::function<Point(Point)> xi_from_tau, Vector Jtau = {}, Vector Jsigma = {});

}
#endif
