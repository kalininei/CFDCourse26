#ifndef CFD_PRINTER_HPP
#define CFD_PRINTER_HPP

#include "cfd/geom/point.hpp"
#include "cfd/mat/i_sparse_mat.hpp"

namespace cfd {

/// @brief debug procedures
namespace dbg {

void ping_printer_cpp();

/**
 * @brief prints sparse matrix to std::cout
 *
 * @param mat  sparse matrix
 *
 * prints '*' for entries that do not present in the stencil.
 */
void print(const ISparseMatrix& mat);

/**
 * @brief prints sparse matrix row to std::cout
 *
 * @param irow  row index
 * @param mat   sparse matrix
 *
 * prints only entries that are in the stencil
 */
void print(size_t irow, const ISparseMatrix& mat);

/**
 * @brief prints sparse matrix row to std::cout
 *
 * @param irow  row index
 * @param mat   sparse matrix
 * @param icol0 start column
 * @param icol1 end column
 *
 * prints only entries that are in the stencil
 */
void print(size_t irow, const ISparseMatrix& mat, size_t col0, size_t col1);

/**
 * @brief prints dense vector to std::cout
 *
 * @param vec  vector to print
 */
void print(const std::vector<double>& vec);

void print(const std::vector<Vector>& vec);

/**
 * @brief prints vector data features
 *
 * @param vec input vector
 */
void print_feat(const std::vector<double>& vec);
} // namespace dbg
} // namespace cfd
#endif
