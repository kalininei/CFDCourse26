#include "lodmat.hpp"

using namespace cfd;

LodMatrix::LodMatrix(size_t n_rows) : data_(n_rows) {}

const std::map<size_t, double>& LodMatrix::row(size_t i_row) const {
    return data_.at(i_row);
}

void LodMatrix::add_value(size_t irow, size_t icol, double value) {
    std::map<size_t, double>& r = data_.at(irow);
    auto found = r.find(icol);
    if (found == r.end()) {
        r.emplace(icol, value);
    } else {
        found->second += value;
    }
}

void LodMatrix::set_value(size_t irow, size_t icol, double value) {
    std::map<size_t, double>& r = data_.at(irow);
    r[icol] = value;
}

void LodMatrix::remove_value(size_t irow, size_t icol) {
    std::map<size_t, double>& r = data_.at(irow);
    auto found = r.find(icol);
    if (found != r.end()) {
        r.erase(found);
    }
}

void LodMatrix::remove_row(size_t irow) {
    data_.at(irow).clear();
}

void LodMatrix::set_unit_row(size_t irow) {
    data_.at(irow) = std::map<size_t, double>{{irow, 1.0}};
}

size_t LodMatrix::n_rows() const {
    return data_.size();
}

double LodMatrix::value(size_t irow, size_t icol) const {
    const std::map<size_t, double>& r = data_.at(irow);
    auto found = r.find(icol);
    if (found == r.end()) {
        return 0;
    } else {
        return found->second;
    }
}

size_t LodMatrix::n_nonzeros() const {
    size_t ret = 0;
    for (const auto& it: data_) {
        ret += it.size();
    }
    return ret;
}

bool LodMatrix::is_in_stencil(size_t irow, size_t icol) const {
    const std::map<size_t, double>& r = data_.at(irow);
    return r.find(icol) != r.end();
}

CsrMatrix LodMatrix::to_csr() const {
    // build csr arrays
    std::vector<size_t> addr{0};
    std::vector<size_t> cols;
    std::vector<double> vals;
    for (size_t irow = 0; irow < n_rows(); ++irow) {
        const std::map<size_t, double>& r = row(irow);
        for (const auto& it: r) {
            cols.push_back(it.first);
            vals.push_back(it.second);
        }
        addr.push_back(cols.size());
    }

    // construct csr matrix
    CsrMatrix ret;
    ret.set_stencil(std::move(addr), std::move(cols));
    ret.set_values(std::move(vals));
    ret.validate();

    return ret;
}

std::vector<double> LodMatrix::mult_vec_p(const double* u) const {
    std::vector<double> ret(n_rows(), 0);

    for (size_t i = 0; i < n_rows(); ++i) {
        for (const auto& it: data_[i]) {
            size_t j = it.first;
            double aij = it.second;
            ret[i] += aij * u[j];
        }
    }

    return ret;
}

double LodMatrix::mult_vec_p(size_t irow, const double* u) const {
    double ret = 0;
    for (const auto& it: data_[irow]) {
        size_t j = it.first;
        double aij = it.second;
        ret += aij * u[j];
    }
    return ret;
}
