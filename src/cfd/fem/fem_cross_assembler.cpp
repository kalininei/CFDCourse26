#include "fem_cross_assembler.hpp"

using namespace cfd;

FemCrossAssembler::FemCrossAssembler(const FemAssembler& fa_rows, const FemAssembler& fa_cols)
    : fa_rows_(fa_rows),
      fa_cols_(fa_cols) {
    // stencil
    std::vector<std::set<size_t>> tab_basis_basis(fa_rows_.n_bases());
    for (size_t ielem = 0; ielem < fa_rows_.n_elements(); ++ielem) {
        for (size_t ibas: fa_rows_.tab_elem_basis(ielem)) {
            for (size_t jbas: fa_cols_.tab_elem_basis(ielem)) {
                tab_basis_basis[ibas].insert(jbas);
            }
        }
    }
    stencil_.set_stencil(tab_basis_basis);

    // element -> csr addresses table
    tab_elem_csr_address_.resize(fa_rows_.n_elements());
    for (size_t ielem = 0; ielem < fa_rows_.n_elements(); ++ielem) {
        size_t n_rows = fa_rows_.element(ielem).basis->size();
        size_t n_cols = fa_cols_.element(ielem).basis->size();
        for (size_t local_row = 0; local_row < n_rows; ++local_row) {
            size_t global_row = fa_rows_.tab_elem_basis(ielem)[local_row];
            for (size_t local_col = 0; local_col < n_cols; ++local_col) {
                size_t global_col = fa_cols_.tab_elem_basis(ielem)[local_col];
                size_t addr = stencil_.get_address(global_row, global_col);
                if (addr == INVALID_INDEX) {
                    _THROW_INTERNAL_ERROR_;
                }
                tab_elem_csr_address_[ielem].push_back(addr);
            }
        }
    }
}

const CsrStencil& FemCrossAssembler::stencil() const {
    return stencil_;
}

void FemCrossAssembler::add_to_global_matrix(double coef, size_t ielem, const std::vector<double>& local_matrix,
                                             std::vector<double>& global_csr_vals) const {
    for (size_t ival = 0; ival < local_matrix.size(); ++ival) {
        size_t a = tab_elem_csr_address_[ielem][ival];
        global_csr_vals[a] += coef * local_matrix[ival];
    }
}
