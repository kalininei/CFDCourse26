#include "i_grid.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Cache
///////////////////////////////////////////////////////////////////////////////
void IGrid::Cache::clear() {
    boundary_faces.clear();
    boundary_points.clear();
    boundary_cells.clear();
    point_cell.clear();
}

void IGrid::Cache::need_boundary_faces(const IGrid& grid) {
    if (boundary_faces.size() > 0) {
        return;
    }
    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        std::array<size_t, 2> cc = grid.tab_face_cell(iface);
        if (cc[0] == INVALID_INDEX || cc[1] == INVALID_INDEX) {
            boundary_faces.push_back(iface);
        }
    }
}

void IGrid::Cache::need_boundary_points(const IGrid& grid) {
    if (boundary_points.size() > 0) {
        return;
    }
    need_boundary_faces(grid);
    std::set<size_t> points;
    for (size_t iface: boundary_faces) {
        for (size_t ipoint: grid.tab_face_point(iface)) {
            points.insert(ipoint);
        }
    }
    boundary_points = std::vector<size_t>(points.begin(), points.end());
}

void IGrid::Cache::need_boundary_cells(const IGrid& grid) {
    if (boundary_cells.size() > 0) {
        return;
    }
    need_boundary_faces(grid);
    std::set<size_t> cells;
    for (size_t iface: boundary_faces) {
        auto [c1, c2] = grid.tab_face_cell(iface);
        if (c1 == INVALID_INDEX) {
            c1 = c2;
        }
        cells.insert(c1);
    }
    boundary_cells = std::vector<size_t>(cells.begin(), cells.end());
}

void IGrid::Cache::need_point_cell(const IGrid& grid) {
    if (point_cell.size() > 0) {
        return;
    }
    point_cell.resize(grid.n_points());
    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        for (size_t ipoint: grid.tab_cell_point(icell)) {
            point_cell[ipoint].push_back(icell);
        }
    }
    for (auto& pc: point_cell) {
        std::sort(pc.begin(), pc.end());
    }
}

///////////////////////////////////////////////////////////////////////////////
// IGrid
///////////////////////////////////////////////////////////////////////////////
std::vector<size_t> IGrid::boundary_faces() const {
    _cache.need_boundary_faces(*this);
    return _cache.boundary_faces;
}

std::vector<size_t> IGrid::boundary_points() const {
    _cache.need_boundary_points(*this);
    return _cache.boundary_points;
}

std::vector<size_t> IGrid::boundary_cells() const {
    _cache.need_boundary_cells(*this);
    return _cache.boundary_cells;
}

std::vector<size_t> IGrid::tab_point_cell(size_t ipoint) const {
    _cache.need_point_cell(*this);
    return _cache.point_cell[ipoint];
}
