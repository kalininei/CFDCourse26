#ifndef __CFD_GEOM_JACOBI_HPP__
#define __CFD_GEOM_JACOBI_HPP__

#include "cfd/cfd_common.hpp"
#include "cfd/geom/point.hpp"

namespace cfd {

struct JacobiMatrix {
    double j11 = 1.0;
    double j12 = 0.0;
    double j13 = 0.0;
    double j21 = 0.0;
    double j22 = 1.0;
    double j23 = 0.0;
    double j31 = 0.0;
    double j32 = 0.0;
    double j33 = 1.0;
    double modj = 1.0;
};

void fill_jacobi_modj(JacobiMatrix &jac);
void fill_jacobi_modj_1d(JacobiMatrix &jac);
void fill_jacobi_modj_2d(JacobiMatrix &jac);

Vector gradient_to_parametric(const JacobiMatrix &jac, Vector grad_x);
Vector gradient_to_parametric_1d(const JacobiMatrix &jac, Vector grad_x);
Vector gradient_to_parametric_2d(const JacobiMatrix &jac, Vector grad_x);

Vector gradient_to_physical(const JacobiMatrix &jac, Vector grad_xi);
Vector gradient_to_physical_1d(const JacobiMatrix &jac, Vector grad_xi);
Vector gradient_to_physical_2d(const JacobiMatrix &jac, Vector grad_xi);

} // namespace cfd
#endif
