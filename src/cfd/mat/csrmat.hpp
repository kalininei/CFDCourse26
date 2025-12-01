#ifndef MAT_CSRMAT_HPP
#define MAT_CSRMAT_HPP

#include "cfd/mat/i_sparse_mat.hpp"

namespace cfd {

class DenseMatrix;

/**
 * @brief 'Compressed sparsed row' stencil
 */
class CsrStencil : public ISparseMatrix {
public:
    CsrStencil() = default;
    CsrStencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols)
        : addr_(std::move(addr)), cols_(std::move(cols)) {}
    virtual ~CsrStencil() = default;
    /**
     * @brief Fills CSR stencil by address and column vectors
     *
     * @param addr  address vector
     * @param cols  column vector
     */
    virtual void set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols);
    virtual void set_stencil(const std::vector<std::set<size_t>>& stencil_set);
    virtual void set_stencil(const CsrStencil& stencil);

    /**
     * @brief Address array getter
     */
    const std::vector<size_t>& addr() const;

    /**
     * @brief Columns array getter
     */
    const std::vector<size_t>& cols() const;

    /**
     * @brief Consistency check. Throws if stencil structure is invalid
     */
    virtual void validate() const;

    size_t get_address(size_t irow, size_t icol) const;

    // overrides
    size_t n_rows() const override;
    size_t n_nonzeros() const override;
    bool is_in_stencil(size_t irow, size_t icol) const override;
    double value(size_t irow, size_t icol) const override;
    std::vector<double> mult_vec_p(const double* u) const override;
    double mult_vec_p(size_t irow, const double* u) const override;

private:
    std::vector<size_t> addr_ = {0};
    std::vector<size_t> cols_;
};

/**
 * @brief 'Compressed sparsed row' matrix format
 */
class CsrMatrix : public CsrStencil {
public:
    CsrMatrix() = default;
    CsrMatrix(const CsrStencil& stencil) : CsrStencil(stencil), vals_(stencil.n_nonzeros(), 0.0) {}
    CsrMatrix(std::vector<size_t>&& addr, std::vector<size_t>&& cols, std::vector<double>&& vals)
        : CsrStencil(std::move(addr), std::move(cols)), vals_(std::move(vals)) {}

    void set_data(std::vector<size_t>&& addr, std::vector<size_t>&& cols, std::vector<double>&& vals);

    void set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols) override;
    void set_stencil(const std::vector<std::set<size_t>>& stencil_set) override;
    void set_stencil(const CsrStencil& stencil) override;

    /**
     * @brief Set matrix values
     *
     * @param vals  values array
     */
    void set_values(std::vector<double>&& vals);

    /**
     * @brief Values array getter
     */
    const std::vector<double>& vals() const;
    std::vector<double>& vals();

    void set_unit_row(size_t irow);
    void set_diagonal(const std::vector<double>& diag_vals);
    double row_sum(size_t irow) const;

    // overrides
    void validate() const override;
    double value(size_t irow, size_t icol) const override;
    std::vector<double> mult_vec_p(const double* u) const override;
    double mult_vec_p(size_t irow, const double* u) const override;

    DenseMatrix to_dense() const;

private:
    std::vector<double> vals_;
};

CsrMatrix assemble_block_matrix(size_t block_n_rows, size_t block_n_cols,
                                const std::vector<std::vector<const CsrMatrix*>>& blocks);

class LodMatrix;
CsrMatrix assemble_block_matrix(size_t block_n_rows, size_t block_n_cols,
                                const std::vector<std::vector<const LodMatrix*>>& blocks);

} // namespace cfd

#endif
