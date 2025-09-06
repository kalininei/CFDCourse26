#ifndef MAT_LODMAT_HPP
#define MAT_LODMAT_HPP

#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/i_sparse_mat.hpp"

namespace cfd {

/**
 * @brief 'List of dictionaries' sparse matrix format
 */
class LodMatrix : public ISparseMatrix {
public:
    /**
     * @brief Construct empty matrix
     *
     * @param n_rows number of rows
     */
    LodMatrix(size_t n_rows);

    /**
     * @brief Get row values
     *
     * @param irow  row index
     * @return column->value dictionary
     */
    const std::map<size_t, double>& row(size_t irow) const;

    /**
     * @brief Adds value to the given matrix entry
     *
     * @param irow   row index
     * @param icol   column index
     * @param value  value to be added
     *
     * Performs Matrix[i, j] += value operation
     */
    void add_value(size_t irow, size_t icol, double value);

    /**
     * @brief Sets value to the given matrix entry
     *
     * @param irow   row index
     * @param icol   column index
     * @param value  value to be set
     *
     * Zero value will not be removed from the stencil.
     * Use LodMatrix::remove_value() function to remove from the stencil
     */
    void set_value(size_t irow, size_t icol, double value);

    /**
     * @brief Removes entry from the matrix stencil
     *
     * @param irow   row index
     * @param icol   column index
     */
    void remove_value(size_t irow, size_t icol);

    /**
     * @brief Removes row from sparse matrix stencil
     *
     * @param irow  row index
     *
     * All values at the given row will be equal to zero
     */
    void remove_row(size_t irow);

    /**
     * @brief sets unit diagonal and zero non-diagonal values for the specified row
     */
    void set_unit_row(size_t irow);

    /**
     * @brief converts to csr matrix format
     *
     * @return csr matrix
     */
    CsrMatrix to_csr() const;

    // overrides
    size_t n_rows() const override;
    double value(size_t irow, size_t icol) const override;
    size_t n_nonzeros() const override;
    bool is_in_stencil(size_t irow, size_t icol) const override;
    std::vector<double> mult_vec_p(const double* u) const override;
    double mult_vec_p(size_t irow, const double* u) const override;

private:
    std::vector<std::map<size_t, double>> data_;
};

} // namespace cfd

#endif
