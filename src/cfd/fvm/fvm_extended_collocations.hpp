#ifndef CFD_FVM_EXTENDED_COLLOCATIONS_HPP
#define CFD_FVM_EXTENDED_COLLOCATIONS_HPP

#include "cfd/grid/i_grid.hpp"

namespace cfd {

///////////////////////////////////////////////////////////////////////////////
// FvmExtendedCollocations
///////////////////////////////////////////////////////////////////////////////

struct FvmExtendedCollocations {
public:
    explicit FvmExtendedCollocations(const IGrid& grid);

    // collocation points
    std::vector<Point> points;

    // indices of cell collocations
    std::vector<size_t> cell_collocations;

    // indices of face collocations
    std::vector<size_t> face_collocations;

    /// number of collocations
    size_t size() const;

    /// index of a cell for the given collocation. Throws if icolloc is not a cell collocation
    size_t cell_index(size_t icolloc) const;

    /// index of a face for the given collocation. Throws if icolloc is not a face collocation
    size_t face_index(size_t icolloc) const;

    // face -> left & right collocation indices
    std::array<size_t, 2> tab_face_colloc(size_t iface) const;

    // gets a collocation index for a boundary face. throws if iface is not boundary
    size_t boundary_colloc(size_t iface) const;

    // colloc -> connected collocs
    std::vector<size_t> tab_colloc_colloc(size_t icolloc) const;

    bool is_boundary_colloc(size_t icolloc) const;
    bool is_internal_colloc(size_t icolloc) const;

private:
    std::vector<std::array<size_t, 2>> tab_face_colloc_;
    std::vector<std::vector<size_t>> tab_colloc_colloc_;
    std::vector<size_t> face_indices_;
};

} // namespace cfd

#endif
