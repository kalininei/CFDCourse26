#ifndef CFD_FVM_DFDN_HPP
#define CFD_FVM_DFDN_HPP

#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/grid/i_grid.hpp"
#include "cfd/mat/lodmat.hpp"

namespace cfd {

///////////////////////////////////////////////////////////////////////////////
// DfDn on faces
///////////////////////////////////////////////////////////////////////////////

struct FvmFacesDn {
    explicit FvmFacesDn(const IGrid& grid);
    FvmFacesDn(const IGrid& grid, const FvmExtendedCollocations& colloc);

    /// computes dfdn for each grid face
    std::vector<double> compute(const std::vector<double>& f) const;
    std::vector<double> compute(const double* f) const;

    /// computes dfdn for the given grid face
    double compute(size_t iface, const std::vector<double>& f) const;
    double compute(size_t iface, const double* f) const;

    /// returns dfdn as a linear combination of collocation values
    const std::map<size_t, double>& linear_combination(size_t iface) const;

private:
    LodMatrix dfdn_;
};

} // namespace cfd
#endif
