#ifndef __CFD_FEM_SORTED_CELL_INFO_HPP__
#define __CFD_FEM_SORTED_CELL_INFO_HPP__

#include "cfd/grid/i_grid.hpp"

namespace cfd {

///////////////////////////////////////////////////////////////////////////////
// 2D Polygon Element Info
///////////////////////////////////////////////////////////////////////////////
struct PolygonElementInfo {
    PolygonElementInfo(const IGrid& grid, size_t icell);

    size_t n_points() const {
        return ipoints.size();
    }

    size_t icell;
    std::vector<size_t> ipoints;
    std::vector<size_t> ifaces;
    std::vector<bool> is_face_reverted;
};

} // namespace cfd

#endif
