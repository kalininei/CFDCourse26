#ifndef MAT_I_SPARSE_MAT_HPP
#define MAT_I_SPARSE_MAT_HPP

#include "cfd/mat/i_mat.hpp"

namespace cfd {

/**
 * @brief Abstract Sparse Matrix interface
 */
class ISparseMatrix : public IMatrix {
public:
    virtual ~ISparseMatrix() = default;

    /**
     * @brief Gets number of entries in the stencil
     */
    virtual size_t n_nonzeros() const = 0;

    /**
     * @brief Checks if entry in non-zero stencil
     *
     * @param irow  row index
     * @param icol  column index
     * @return true if [irow, icol] is in nonzero stencil
     */
    virtual bool is_in_stencil(size_t irow, size_t icol) const = 0;
};

} // namespace cfd

#endif
