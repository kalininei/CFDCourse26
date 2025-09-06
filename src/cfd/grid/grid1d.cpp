#include "grid1d.hpp"
#include "cfd/grid/vtk.hpp"

using namespace cfd;

Grid1D::Grid1D(double left, double right, size_t n_cells) {
    points_.push_back(Point(left));
    double h = (right - left) / n_cells;
    for (size_t i = 0; i < n_cells; ++i) {
        points_.push_back(Point(points_.back()[0] + h));
    }
}

size_t Grid1D::n_points() const {
    return points_.size();
}

size_t Grid1D::n_cells() const {
    return points_.size() - 1;
}

size_t Grid1D::n_faces() const {
    return points_.size();
}

Point Grid1D::cell_center(size_t icell) const {
    return (points_.at(icell) + points_.at(icell + 1)) / 2.0;
}

double Grid1D::cell_volume(size_t icell) const {
    return points_.at(icell + 1).x() - points_.at(icell).x();
}

Vector Grid1D::face_normal(size_t iface) const {
    return Vector(1.0, 0.0, 0.0);
}

double Grid1D::face_area(size_t iface) const {
    return 1.0;
}

Point Grid1D::face_center(size_t iface) const {
    return points_.at(iface);
}

Point Grid1D::point(size_t ipoint) const {
    return points_.at(ipoint);
}

std::vector<Point> Grid1D::points() const {
    return points_;
}

std::vector<size_t> Grid1D::tab_cell_point(size_t icell) const {
    return {icell, icell + 1};
}

std::vector<size_t> Grid1D::tab_face_point(size_t iface) const {
    return {iface};
}

std::array<size_t, 2> Grid1D::tab_face_cell(size_t iface) const {
    if (iface == 0) {
        return {INVALID_INDEX, 0};
    } else if (iface == points_.size() - 1) {
        return {iface - 1, INVALID_INDEX};
    } else {
        return {iface - 1, iface};
    }
}

std::vector<size_t> Grid1D::tab_cell_face(size_t icell) const {
    return {icell, icell + 1};
}

void Grid1D::save_vtk(std::string fname) const {
    std::ofstream fs(fname);
    VtkUtils::append_header("Grid1", fs);
    VtkUtils::append_points(this->points(), fs);
    // Cells
    fs << "CELLS  " << this->n_cells() << "   " << 3 * this->n_cells() << std::endl;
    for (size_t i = 0; i < this->n_cells(); ++i) {
        auto v = this->tab_cell_point(i);
        fs << 2 << " " << v[0] << " " << v[1] << std::endl;
    }
    fs << "CELL_TYPES  " << this->n_cells() << std::endl;
    for (size_t i = 0; i < this->n_cells(); ++i)
        fs << 3 << std::endl;
}
