#ifndef CFD_DENSE_MAT_HPP
#define CFD_DENSE_MAT_HPP

#include "cfd/mat/i_mat.hpp"

namespace cfd {

class CsrMatrix;

class DenseMatrix : public IMatrix {
public:
    DenseMatrix(size_t nrows, size_t ncols);
    DenseMatrix(size_t nrows, size_t ncols, const std::vector<double>& values);

    void set_value(size_t irow, size_t icol, double value);

    DenseMatrix transpose() const;
    DenseMatrix mult_mat(const DenseMatrix& mat) const;
    DenseMatrix inverse() const;

    size_t n_cols() const;
    const std::vector<double>& vals() const;

    // overriden
    size_t n_rows() const override;
    double value(size_t irow, size_t icol) const override;
    std::vector<double> mult_vec_p(const double* u) const override;
    double mult_vec_p(size_t irow, const double* u) const override;

    CsrMatrix to_csr() const;

private:
    const size_t nrows_;
    const size_t ncols_;
    std::vector<double> data_;
};

} // namespace cfd
#endif
