#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP

#include "cfd/mat/csrmat.hpp"

namespace cfd {

/**
 * @brief Amgcl Solver for csr matrices
 */
class AmgcMatrixSolver {
public:
    /**
     * @param maxit  maximum allowed number of iterations
     * @param eps    tolerance value
     */
    AmgcMatrixSolver(int maxit = 1000, double eps = 1e-8);

    /**
     * @param amgc_params  set of AMGCL non-default parameters
     */
    AmgcMatrixSolver(std::initializer_list<std::pair<std::string, std::string>> amgc_params);

    ~AmgcMatrixSolver();

    /**
     * @brief Sets target matrix
     *
     * @param mat  target matrix
     *
     * Matrix will be copied to the internal structure and can be destroyed
     */
    void set_matrix(const CsrMatrix& mat);

    /**
     * @brief Sets target matrix
     *
     * @param mat_stencil  csr matrix stencil
     * @param mat_values   matrix values
     *
     * Matrix will be copied to the internal structure and can be destroyed
     */
    void set_matrix(const CsrStencil& mat_stencil, const std::vector<double>& mat_values);

    /**
     * @brief Solves slae Ax = rhs
     *
     * @param      rhs  right hand side vector
     * @param[out] ret  output vector
     *
     * Given values of output vector will be used as the initial values.
     * If it is empty, solution will be initialized with zeros.
     */
    void solve(const std::vector<double>& rhs, std::vector<double>& ret) const;

    static void solve_slae(const CsrMatrix& mat, const std::vector<double>& rhs, std::vector<double>& x,
                           int maxit = 1000, double eps = 1e-8);

private:
    int maxit_;
    double tolerance_;
    std::map<std::string, std::string> params_;

    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace cfd

#endif
