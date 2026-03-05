#include "grid_manip.hpp"
#include "cfd/geom/point.hpp"

using namespace cfd;

UnstructuredGrid2D dbg::to_unstructured(const IGrid2D& grid) {
    std::vector<std::vector<size_t>> cp_tab;
    for (size_t i = 0; i < grid.n_cells(); ++i) {
        cp_tab.push_back(grid.tab_cell_point(i));
    }
    return UnstructuredGrid2D(grid.points(), cp_tab);
}

UnstructuredGrid2D dbg::rotate_grid(const IGrid2D& grid, double angle_radians, Point center) {
    std::vector<std::vector<size_t>> cp_tab;
    for (size_t i = 0; i < grid.n_cells(); ++i) {
        cp_tab.push_back(grid.tab_cell_point(i));
    }
    UnstructuredGrid2D g1(grid.points(), cp_tab);

    return g1.copy_modify([center, angle_radians](Point p) { return rotate(p, angle_radians, center); });
}
