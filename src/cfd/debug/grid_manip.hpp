#ifndef CFD_COURSE_GRID_MANIP_HPP
#define CFD_COURSE_GRID_MANIP_HPP

#include "cfd/grid/unstructured_grid2d.hpp"

namespace cfd {
namespace dbg {

UnstructuredGrid2D to_unstructured(const IGrid2D& grid);
UnstructuredGrid2D rotate_grid(const IGrid2D& grid, double angle_radians, Point center = {});

} // namespace dbg
} // namespace cfd

#endif
