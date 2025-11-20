#include "csrmat.hpp"
#include "lodmat.hpp"
#include <numeric>

using namespace cfd;

void CsrStencil::set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols) {
    addr_ = std::move(addr);
    cols_ = std::move(cols);
}

void CsrStencil::set_stencil(const std::vector<std::set<size_t>>& stencil_set) {
    addr_ = std::vector<size_t>(1, 0);
    cols_.clear();

    for (size_t irow = 0; irow < stencil_set.size(); ++irow) {
        const std::set<size_t>& cols = stencil_set[irow];
        addr_.push_back(addr_.back() + cols.size());
        for (size_t col: cols) {
            cols_.push_back(col);
        }
    }
}

void CsrStencil::set_stencil(const CsrStencil& stencil) {
    std::vector<size_t> addr = stencil.addr();
    std::vector<size_t> cols = stencil.cols();
    set_stencil(std::move(addr), std::move(cols));
}

size_t CsrStencil::n_nonzeros() const {
    return cols_.size();
}

size_t CsrStencil::n_rows() const {
    return addr_.size() - 1;
}

const std::vector<size_t>& CsrStencil::addr() const {
    return addr_;
}

const std::vector<size_t>& CsrStencil::cols() const {
    return cols_;
}

void CsrStencil::validate() const {
    // sizes
    if (addr_.size() < 1) {
        throw std::runtime_error("addr array should have more then zero entries");
    }
    if (cols_.size() != addr_.back()) {
        throw std::runtime_error("cols array size should match last addr entry");
    }
    // non-decreasing
    for (size_t i = 1; i < addr_.size(); ++i) {
        if (addr_[i] < addr_[i - 1]) {
            throw std::runtime_error("addr array should be non-decreasing");
        }
    }
}

bool CsrStencil::is_in_stencil(size_t irow, size_t icol) const {
    size_t start = addr_.at(irow);
    size_t end = addr_.at(irow + 1);
    for (size_t i = start; i < end; ++i) {
        if (cols_.at(i) == icol) {
            return true;
        }
    }
    return false;
}

double CsrStencil::value(size_t, size_t) const {
    throw std::runtime_error("CsrStencil has no values");
}

std::vector<double> CsrStencil::mult_vec_p(const double*) const {
    throw std::runtime_error("CsrStencil has no values");
}

double CsrStencil::mult_vec_p(size_t, const double*) const {
    throw std::runtime_error("CsrStencil has no values");
}

size_t CsrStencil::get_address(size_t irow, size_t icol) const {
    std::vector<size_t>::const_iterator it_start = cols_.begin() + addr_.at(irow);
    std::vector<size_t>::const_iterator it_end = cols_.begin() + addr_.at(irow + 1);
    auto fnd = std::lower_bound(it_start, it_end, icol);
    if (fnd != it_end && *fnd == icol) {
        size_t a = fnd - cols_.begin();
        return a;
    }
    return INVALID_INDEX;
}

void CsrMatrix::set_values(std::vector<double>&& vals) {
    vals_ = std::move(vals);
}

const std::vector<double>& CsrMatrix::vals() const {
    return vals_;
}

std::vector<double>& CsrMatrix::vals() {
    return vals_;
}

void CsrMatrix::validate() const {
    CsrStencil::validate();

    // values size
    if (vals_.size() != n_nonzeros()) {
        throw std::runtime_error("values array should have same size as the columns arrays");
    }
}

double CsrMatrix::value(size_t irow, size_t icol) const {
    size_t a = get_address(irow, icol);
    if (a != INVALID_INDEX) {
        return vals_[a];
    } else {
        return 0.0;
    }
}

std::vector<double> CsrMatrix::mult_vec_p(const double* u) const {
    const std::vector<size_t>& a = addr();
    const std::vector<size_t>& c = cols();
    const std::vector<double>& v = vals();

    std::vector<double> ret(n_rows(), 0);
    for (size_t irow = 0; irow < n_rows(); ++irow) {
        size_t start = a[irow];
        size_t end = a[irow + 1];
        for (size_t i = start; i < end; ++i) {
            ret[irow] += v[i] * u[c[i]];
        }
    }

    return ret;
}

double CsrMatrix::mult_vec_p(size_t irow, const double* u) const {
    const std::vector<size_t>& a = addr();
    const std::vector<size_t>& c = cols();
    const std::vector<double>& v = vals();

    double ret = 0;
    size_t start = a.at(irow);
    size_t end = a.at(irow + 1);
    for (size_t i = start; i < end; ++i) {
        ret += v[i] * u[c[i]];
    }
    return ret;
}

void CsrMatrix::set_unit_row(size_t irow) {
    const std::vector<size_t>& a = addr();
    const std::vector<size_t>& c = cols();

    const size_t start = a.at(irow);
    const size_t end = a.at(irow + 1);
    for (size_t i = start; i < end; ++i) {
        vals_[i] = (c[i] == irow) ? 1.0 : 0.0;
    }
}

void CsrMatrix::set_diagonal(const std::vector<double>& diag_vals) {
    const std::vector<size_t>& a = addr();
    const std::vector<size_t>& c = cols();

    for (size_t irow = 0; irow < n_rows(); ++irow) {
        const size_t start = a.at(irow);
        const size_t end = a.at(irow + 1);
        for (size_t i = start; i < end; ++i) {
            if (c[i] == irow) {
                vals_[i] = diag_vals[irow];
            }
        }
    }
}

CsrMatrix cfd::assemble_block_matrix(size_t block_n_rows, size_t block_n_cols,
                                     const std::vector<std::vector<const CsrMatrix*>>& blocks) {

    size_t nrows = blocks.size() * block_n_rows;
    std::vector<size_t> cols;
    std::vector<double> vals;
    std::vector<size_t> n_row_nonzeros(nrows, 0);

    size_t col_margin = 0;
    size_t row_margin = 0;
    for (size_t i_block_row = 0; i_block_row < blocks.size(); ++i_block_row) {

        for (size_t irow = 0; irow < block_n_rows; ++irow) {

            for (size_t i_block_col = 0; i_block_col < blocks[i_block_row].size(); ++i_block_col) {
                const CsrMatrix* block = blocks[i_block_row][i_block_col];

                if (block) {
                    // n_nonzeros in row
                    size_t nz = block->addr()[irow + 1] - block->addr()[irow];
                    n_row_nonzeros[irow + row_margin] += nz;

                    // values
                    const double* v = &block->vals()[block->addr()[irow]];
                    vals.insert(vals.end(), v, v + nz);

                    // columns
                    const size_t* c = &block->cols()[block->addr()[irow]];
                    for (size_t i = 0; i < nz; ++i)
                        cols.push_back(c[i] + col_margin);
                }

                col_margin += block_n_cols;
            }
            col_margin = 0;
        }
        row_margin += block_n_rows;
    }

    // assemble addr with running sum
    std::vector<size_t> addr(nrows + 1, 0);
    std::partial_sum(n_row_nonzeros.begin(), n_row_nonzeros.end(), addr.begin() + 1);

    return CsrMatrix(std::move(addr), std::move(cols), std::move(vals));
}

CsrMatrix cfd::assemble_block_matrix(size_t n_block_rows, size_t n_block_cols,
                                     const std::vector<std::vector<const LodMatrix*>>& blocks) {
    std::vector<CsrMatrix> data;
    size_t nblocks = 0;
    for (const auto& r: blocks) {
        nblocks += r.size();
    };
    data.reserve(nblocks);

    std::vector<std::vector<const CsrMatrix*>> ret(blocks.size());

    for (size_t i = 0; i < blocks.size(); ++i) {
        for (size_t j = 0; j < blocks[i].size(); ++j) {
            const LodMatrix* block = blocks[i][j];

            if (block) {
                data.push_back(block->to_csr());
                ret[i].push_back(&data.back());
            } else {
                ret[i].push_back(nullptr);
            }
        }
    }

    return assemble_block_matrix(n_block_rows, n_block_cols, ret);
}

void CsrMatrix::set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols) {
    CsrStencil::set_stencil(std::move(addr), std::move(cols));
    vals_ = std::vector<double>(CsrStencil::cols().size(), 0.0);
}

void CsrMatrix::set_stencil(const std::vector<std::set<size_t>>& stencil_set) {
    CsrStencil::set_stencil(stencil_set);
    vals_ = std::vector<double>(CsrStencil::cols().size(), 0.0);
}

void CsrMatrix::set_stencil(const CsrStencil& stencil) {
    CsrStencil::set_stencil(stencil);
    vals_ = std::vector<double>(CsrStencil::cols().size(), 0.0);
}

double CsrMatrix::row_sum(size_t irow) const {
    double ret = 0;
    const std::vector<size_t>& a = addr();

    const size_t start = a.at(irow);
    const size_t end = a.at(irow + 1);
    for (size_t i = start; i < end; ++i) {
        ret += vals_[i];
    }

    return ret;
}
