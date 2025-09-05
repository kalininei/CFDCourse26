#ifndef __CFD_VECMAT_HPP__
#define __CFD_VECMAT_HPP__

#include "cfd/geom/point.hpp"
#include "cfd/mat/csrmat.hpp"
#include <vector>

// -> abs(m*u - rhs)
std::vector<double> compute_residual_vec(const cfd::CsrMatrix &m, const std::vector<double> &rhs, const std::vector<double> &u);
// -> max(abs(m*u - rhs))
double compute_residual(const cfd::CsrMatrix &m, const std::vector<double> &rhs, const std::vector<double> &u);

// -> v1 + coef*v2
std::vector<double> vector_sum(
    const std::vector<double> &v1,
    double coef,
    const std::vector<double> &v2);

std::vector<cfd::Vector> vector_sum(
    const std::vector<cfd::Vector> &v1,
    double coef,
    const std::vector<cfd::Vector> &v2);

double harmonic_mean(double a, double b);

#endif
