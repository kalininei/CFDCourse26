#ifndef CFD_FEM_CROSS_ASSEMBLER_HPP
#define CFD_FEM_CROSS_ASSEMBLER_HPP

#include "cfd/fem/fem_assembler.hpp"

namespace cfd{

class FemCrossAssembler{
public:
    FemCrossAssembler(const FemAssembler& fa_rows, const FemAssembler& fa_cols);
    const CsrStencil& stencil() const;
    void add_to_global_matrix(double coef, size_t ielem, const std::vector<double>& local_matrix,
                              std::vector<double>& global_csr_vals) const;
private:
    const FemAssembler& fa_rows_;
    const FemAssembler& fa_cols_;

    CsrStencil stencil_;
    std::vector<std::vector<size_t>> tab_elem_csr_address_;
};


}
#endif
