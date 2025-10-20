#include "fem_sorted_cell_info.hpp"

using namespace cfd;

PolygonElementInfo::PolygonElementInfo(const IGrid& grid, size_t icell)
    : icell(icell), ipoints(grid.tab_cell_point(icell)) {
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
}
