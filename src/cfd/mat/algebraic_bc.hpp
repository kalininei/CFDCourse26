#ifndef __CFD_MAT_ALGEBRAIC_BC_HPP__
#define __CFD_MAT_ALGEBRAIC_BC_HPP__

#include "cfd/mat/csrmat.hpp"

namespace cfd {

// Dirichlet
// vals: [i -> f(x_i)] dirichlet values
void algebraic_bc_dirichlet(const std::map<size_t, double>& vals, CsrMatrix& lhs);
void algebraic_bc_dirichlet(const std::map<size_t, double>& vals, std::vector<double>& rhs);

// Periodic
// connections: [i -> j] periodic pairs:
void algebraic_bc_periodic(const std::map<size_t, size_t>& connections, CsrMatrix& lhs);
void algebraic_bc_periodic(const std::map<size_t, size_t>& connections, std::vector<double>& rhs);

}; // namespace cfd

#endif
