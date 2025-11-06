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

    // if this is a boundary cell, reverts ipoints, ifaces order so that the boundary face come first
    // returns true if renumbering was done
    bool start_from_boundary();

    size_t icell;
    std::vector<size_t> ipoints;
    std::vector<size_t> ifaces;
    std::vector<bool> is_face_reverted;
    bool is_boundary_cell;
    std::vector<bool> is_boundary_face;
};

} // namespace cfd

#endif
