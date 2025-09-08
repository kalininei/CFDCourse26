#include "densemat.hpp"

using namespace cfd;

DenseMatrix::DenseMatrix(size_t nrows, size_t ncols) : nrows_(nrows), ncols_(ncols), data_(nrows * ncols, 0) {}

DenseMatrix::DenseMatrix(size_t nrows, size_t ncols, const std::vector<double>& values)
    : nrows_(nrows), ncols_(ncols), data_(values) {}

void DenseMatrix::set_value(size_t irow, size_t icol, double value) {
    data_[irow * ncols_ + icol] = value;
}

DenseMatrix DenseMatrix::transpose() const {
    DenseMatrix ret(ncols_, nrows_);

    for (size_t i = 0; i < nrows_; ++i)
        for (size_t j = 0; j < ncols_; ++j) {
            size_t old_index = i * ncols_ + j;
            size_t new_index = j * nrows_ + i;

            ret.data_[new_index] = data_[old_index];
        }

    return ret;
}

DenseMatrix DenseMatrix::mult_mat(const DenseMatrix& mat) const {
    if (n_cols() != mat.n_rows()) {
        _THROW_INTERNAL_ERROR_;
    }
    DenseMatrix ret(n_rows(), mat.n_cols());
    for (size_t i = 0; i < n_rows(); ++i)
        for (size_t j = 0; j < mat.n_cols(); ++j) {
            double sum = 0;
            for (size_t k = 0; k < n_cols(); ++k) {
                sum += value(i, k) * mat.value(k, j);
            }
            ret.set_value(i, j, sum);
        }

    return ret;
}

DenseMatrix DenseMatrix::inverse() const {
    if (nrows_ == 1) {
        return DenseMatrix(1, 1, {1.0 / data_[0]});
    } else if (nrows_ == 2) {
        double det = data_[0] * data_[3] - data_[1] * data_[2];
        return DenseMatrix(2, 2, {data_[3] / det, -data_[2] / det, -data_[1] / det, data_[0] / det});
    } else if (nrows_ == 3) {
        _THROW_NOT_IMP_;
    } else {
        _THROW_NOT_IMP_;
    }
}

size_t DenseMatrix::n_cols() const {
    return ncols_;
}

const std::vector<double>& DenseMatrix::vals() const {
    return data_;
}

size_t DenseMatrix::n_rows() const {
    return nrows_;
}

double DenseMatrix::value(size_t irow, size_t icol) const {
    return data_[irow * ncols_ + icol];
}

std::vector<double> DenseMatrix::mult_vec_p(const double*) const {
    _THROW_NOT_IMP_;
}

double DenseMatrix::mult_vec_p(size_t, const double*) const {
    _THROW_NOT_IMP_;
}
