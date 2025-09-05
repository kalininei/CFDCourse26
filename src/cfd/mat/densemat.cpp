#include "densemat.hpp"

using namespace cfd;

DenseMatrix::DenseMatrix(size_t nrows, size_t ncols)
    : _nrows(nrows), _ncols(ncols), _data(nrows * ncols, 0) {}

DenseMatrix::DenseMatrix(size_t nrows, size_t ncols, const std::vector<double> &values)
    : _nrows(nrows), _ncols(ncols), _data(values) {}

void DenseMatrix::set_value(size_t irow, size_t icol, double value) {
    _data[irow * _ncols + icol] = value;
}

DenseMatrix DenseMatrix::transpose() const {
    DenseMatrix ret(_ncols, _nrows);

    for (size_t i = 0; i < _nrows; ++i)
        for (size_t j = 0; j < _ncols; ++j) {
            size_t old_index = i * _ncols + j;
            size_t new_index = j * _nrows + i;

            ret._data[new_index] = _data[old_index];
        }

    return ret;
}

DenseMatrix DenseMatrix::mult_mat(const DenseMatrix &mat) const {
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
    if (_nrows == 1) {
        return DenseMatrix(1, 1, {1.0 / _data[0]});
    } else if (_nrows == 2) {
        double det = _data[0] * _data[3] - _data[1] * _data[2];
        return DenseMatrix(2, 2, {_data[3] / det, -_data[2] / det, -_data[1] / det, _data[0] / det});
    } else if (_nrows == 3) {
        _THROW_NOT_IMP_;
    } else {
        _THROW_NOT_IMP_;
    }
}

size_t DenseMatrix::n_cols() const {
    return _ncols;
}

const std::vector<double> &DenseMatrix::vals() const {
    return _data;
}

size_t DenseMatrix::n_rows() const {
    return _nrows;
}

double DenseMatrix::value(size_t irow, size_t icol) const {
    return _data[irow * _ncols + icol];
}

std::vector<double> DenseMatrix::mult_vec_p(const double *u) const {
    _THROW_NOT_IMP_;
}

double DenseMatrix::mult_vec_p(size_t irow, const double *u) const {
    _THROW_NOT_IMP_;
}
