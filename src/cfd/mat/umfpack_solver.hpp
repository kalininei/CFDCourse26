#ifndef CFD_SUPERLU_SOLVER_HPP
#define CFD_SUPERLU_SOLVER_HPP

#include "cfd/mat/csrmat.hpp"

namespace cfd{

class UmfpackSolver {
public:
    UmfpackSolver();
    ~UmfpackSolver();
    
    void set_matrix(const CsrMatrix& m);
    void solve(const std::vector<double>& rhs, std::vector<double>& x);

    static void solve_slae(const CsrMatrix& mat, const std::vector<double>& rhs, std::vector<double>& x);
private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};



}
#endif
