#include "fem_sorted_cell_info.hpp"
#include <algorithm>

using namespace cfd;

PolygonElementInfo::PolygonElementInfo(const IGrid& grid, size_t icell1)
    : icell(icell1), ipoints(grid.tab_cell_point(icell)) {
    std::vector<size_t> orig_ifaces = grid.tab_cell_face(icell);

    for (size_t ip = 0; ip < ipoints.size(); ++ip) {
        size_t p0 = ipoints[ip];
        size_t p1 = ipoints[(ip + 1) % ipoints.size()];

        for (size_t iface : orig_ifaces) {
            std::vector<size_t> face_points = grid.tab_face_point(iface);
            if (face_points.size() != 2) {
                _THROW_INTERNAL_ERROR_;
            }
            if (p0 == face_points[0] && p1 == face_points[1]) {
                ifaces.push_back(iface);
                is_face_reverted.push_back(false);
                break;
            } else if (p0 == face_points[1] && p1 == face_points[0]) {
                ifaces.push_back(iface);
                is_face_reverted.push_back(true);
                break;
            }
        }
    }

    if (ipoints.size() != ifaces.size()) {
        _THROW_INTERNAL_ERROR_;
    }

    for (size_t iface : ifaces) {
        auto [c1, c2] = grid.tab_face_cell(iface);
        is_boundary_face.push_back(c1 == INVALID_INDEX || c2 == INVALID_INDEX);
    }
    is_boundary_cell = std::any_of(is_boundary_face.begin(), is_boundary_face.end(), [](bool f) { return f; });
}

bool PolygonElementInfo::start_from_boundary() {
    auto fnd = std::find(is_boundary_face.begin(), is_boundary_face.end(), true);
    if (fnd != is_boundary_face.begin() && fnd != is_boundary_face.end()) {
        size_t mid = fnd - is_boundary_face.begin();
        std::rotate(ipoints.begin(), ipoints.begin() + mid, ipoints.end());
        std::rotate(ifaces.begin(), ifaces.begin() + mid, ifaces.end());
        std::rotate(is_boundary_face.begin(), is_boundary_face.begin() + mid, is_boundary_face.end());
        return true;
    }
    return false;
}
