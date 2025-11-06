#include "fvm_extended_collocations.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// FvmExtendedCollocations
///////////////////////////////////////////////////////////////////////////////

FvmExtendedCollocations::FvmExtendedCollocations(const IGrid& grid) {
    // cell loop
    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        // collocations at the cell center
        points.push_back(grid.cell_center(icell));
        cell_collocations.push_back(icell);
    }

    // face loop
    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        std::array<size_t, 2> cells = grid.tab_face_cell(iface);

        // face -> collocation-points connectivity
        tab_face_colloc_.push_back({cells[0], cells[1]});

        // collocations at the boundary faces
        if (cells[0] == INVALID_INDEX || cells[1] == INVALID_INDEX) {
            points.push_back(grid.face_center(iface));
            face_collocations.push_back(points.size() - 1);
            face_indices_.push_back(iface);

            if (cells[0] == INVALID_INDEX) {
                tab_face_colloc_.back()[0] = points.size() - 1;
            } else {
                tab_face_colloc_.back()[1] = points.size() - 1;
            }
        }
    }

    // collocations connectivity
    tab_colloc_colloc_.resize(points.size());
    for (const std::array<size_t, 2>& fc: tab_face_colloc_) {
        tab_colloc_colloc_[fc[1]].push_back(fc[0]);
        tab_colloc_colloc_[fc[0]].push_back(fc[1]);
    }
}

size_t FvmExtendedCollocations::size() const {
    return points.size();
}

size_t FvmExtendedCollocations::face_index(size_t icolloc) const {
    if (icolloc < cell_collocations.size()) {
        _THROW_INTERNAL_ERROR_;
    } else {
        return face_indices_[icolloc - cell_collocations.size()];
    }
}

size_t FvmExtendedCollocations::cell_index(size_t icolloc) const {
    if (icolloc >= cell_collocations.size()) {
        _THROW_INTERNAL_ERROR_;
    } else {
        return icolloc;
    }
}

std::array<size_t, 2> FvmExtendedCollocations::tab_face_colloc(size_t iface) const {
    return tab_face_colloc_[iface];
}

std::vector<size_t> FvmExtendedCollocations::tab_colloc_colloc(size_t icolloc) const {
    return tab_colloc_colloc_[icolloc];
}

bool FvmExtendedCollocations::is_boundary_colloc(size_t icolloc) const {
    return icolloc >= cell_collocations.size();
}

bool FvmExtendedCollocations::is_internal_colloc(size_t icolloc) const {
    return icolloc < cell_collocations.size();
}

size_t FvmExtendedCollocations::boundary_colloc(size_t iface) const {
    std::array<size_t, 2> ij = tab_face_colloc(iface);
    if (is_boundary_colloc(ij[0])) {
        return ij[0];
    } else if (is_boundary_colloc(ij[1])) {
        return ij[1];
    } else {
        throw std::runtime_error("not a boundary face: " + std::to_string(iface));
    }
}
