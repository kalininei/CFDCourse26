#include "regular_grid2d.hpp"
#include "cfd/grid/vtk.hpp"

using namespace cfd;

RegularGrid2D::RegularGrid2D(double x0, double x1, double y0, double y1, size_t nx, size_t ny) {
    x_.push_back(x0);
    double hx = (x1 - x0) / (double)nx;
    for (size_t i = 0; i < nx; ++i) {
        x_.push_back(x_.back() + hx);
    }
    y_.push_back(y0);
    double hy = (y1 - y0) / (double)ny;
    for (size_t i = 0; i < ny; ++i) {
        y_.push_back(y_.back() + hy);
    }
    actnum_.resize(n_cells(), 1);
    set_face_types();
}

RegularGrid2D::RegularGrid2D(const std::vector<double>& x, const std::vector<double>& y) : x_(x), y_(y) {

    actnum_.resize(n_cells(), 1);
}

double RegularGrid2D::Lx() const {
    return x_.back() - x_[0];
}

double RegularGrid2D::Ly() const {
    return y_.back() - y_[0];
}

size_t RegularGrid2D::nx() const {
    return x_.size() - 1;
}

size_t RegularGrid2D::ny() const {
    return y_.size() - 1;
}

size_t RegularGrid2D::n_points() const {
    return x_.size() * y_.size();
}

size_t RegularGrid2D::n_cells() const {
    return (x_.size() - 1) * (y_.size() - 1);
}

size_t RegularGrid2D::n_faces() const {
    return (nx() + 1) * ny() + nx() * (ny() + 1);
}

Point RegularGrid2D::point(size_t ipoint) const {
    size_t ix = ipoint % x_.size();
    size_t iy = ipoint / x_.size();
    return Point(x_[ix], y_[iy]);
}

Point RegularGrid2D::cell_center(size_t icell) const {
    size_t ix = icell % (x_.size() - 1);
    size_t iy = icell / (x_.size() - 1);
    return Point(0.5 * (x_[ix] + x_[ix + 1]), 0.5 * (y_[iy] + y_[iy + 1]));
}

double RegularGrid2D::cell_volume(size_t icell) const {
    size_t ix = icell % (x_.size() - 1);
    size_t iy = icell / (x_.size() - 1);
    return (x_[ix + 1] - x_[ix]) * (y_[iy + 1] - y_[iy]);
}

Vector RegularGrid2D::face_normal(size_t iface) const {
    size_t n_xfaces = (ny() + 1) * nx();
    if (iface < n_xfaces) {
        return Vector(0.0, 1.0);
    } else {
        return Vector(1.0, 0.0);
    }
}

double RegularGrid2D::face_area(size_t iface) const {
    size_t n_xfaces = (ny() + 1) * nx();
    if (iface < n_xfaces) {
        size_t ix = iface % nx();
        return x_[ix + 1] - x_[ix];
    } else {
        size_t iy = (iface - n_xfaces) / (nx() + 1);
        return y_[iy + 1] - y_[iy];
    }
}

Point RegularGrid2D::face_center(size_t iface) const {
    size_t n_xfaces = (ny() + 1) * nx();
    if (iface < n_xfaces) {
        size_t ix = iface % nx();
        size_t iy = iface / nx();
        return {0.5 * (x_[ix + 1] + x_[ix]), y_[iy]};
    } else {
        size_t ix = (iface - n_xfaces) % (nx() + 1);
        size_t iy = (iface - n_xfaces) / (nx() + 1);
        return {x_[ix], 0.5 * (y_[iy + 1] + y_[iy])};
    }
}

std::vector<Point> RegularGrid2D::points() const {
    std::vector<Point> ret;
    for (size_t i = 0; i < n_points(); ++i) {
        ret.push_back(point(i));
    }
    return ret;
}

std::vector<size_t> RegularGrid2D::tab_cell_point(size_t icell) const {
    size_t ix = icell % (x_.size() - 1);
    size_t iy = icell / (x_.size() - 1);

    size_t n0 = ix + iy * x_.size();
    size_t n1 = ix + 1 + iy * x_.size();
    size_t n2 = ix + 1 + (iy + 1) * x_.size();
    size_t n3 = ix + (iy + 1) * x_.size();
    return {n0, n1, n2, n3};
}

std::array<size_t, 2> RegularGrid2D::tab_face_cell(size_t iface) const {
    size_t n_xfaces = (ny() + 1) * nx();
    if (iface < n_xfaces) {
        size_t ix = iface % nx();
        size_t iy = iface / nx();
        if (iy == 0) {
            return {INVALID_INDEX, iy * nx() + ix};
        } else if (iy == ny()) {
            return {(iy - 1) * nx() + ix, INVALID_INDEX};
        } else {
            return {(iy - 1) * nx() + ix, iy * nx() + ix};
        }
    } else {
        size_t ix = (iface - n_xfaces) % (nx() + 1);
        size_t iy = (iface - n_xfaces) / (nx() + 1);
        if (ix == 0) {
            return {INVALID_INDEX, iy * nx() + ix};
        } else if (ix == nx()) {
            return {iy * nx() + ix - 1, INVALID_INDEX};
        } else {
            return {iy * nx() + ix - 1, iy * nx() + ix};
        }
    }
}

std::vector<size_t> RegularGrid2D::tab_face_point(size_t iface) const {
    size_t n_xfaces = (ny() + 1) * nx();
    if (iface < n_xfaces) {
        size_t ix = iface % nx();
        size_t iy = iface / nx();
        return {iy * (nx() + 1) + ix + 1, iy * (nx() + 1) + ix};
    } else {
        size_t ix = (iface - n_xfaces) % (nx() + 1);
        size_t iy = (iface - n_xfaces) / (nx() + 1);
        return {iy * (nx() + 1) + ix, (iy + 1) * (nx() + 1) + ix};
    }
}

std::vector<size_t> RegularGrid2D::tab_cell_face(size_t icell) const {
    size_t ix = icell % (x_.size() - 1);
    size_t iy = icell / (x_.size() - 1);
    size_t n_xfaces = (ny() + 1) * nx();

    return {iy * nx() + ix, n_xfaces + iy * (nx() + 1) + ix + 1, (iy + 1) * nx() + ix, n_xfaces + iy * (nx() + 1) + ix};
}

void RegularGrid2D::save_vtk(std::string fname) const {
    std::ofstream fs(fname);
    VtkUtils::append_header("Grid2", fs);
    VtkUtils::append_points(this->points(), fs);
    // Cells
    fs << "CELLS  " << this->n_cells() << "   " << 5 * this->n_cells() << std::endl;
    for (size_t i = 0; i < this->n_cells(); ++i) {
        auto v = this->tab_cell_point(i);
        fs << 4 << " " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;
    }
    fs << "CELL_TYPES  " << this->n_cells() << std::endl;
    for (size_t i = 0; i < this->n_cells(); ++i)
        fs << 9 << std::endl;
}

size_t RegularGrid2D::to_linear_point_index(split_index_t point_split_index) const {
    return point_split_index[0] + x_.size() * point_split_index[1];
}

auto RegularGrid2D::to_split_point_index(size_t ipoint) const -> split_index_t {
    return split_index_t{ipoint % x_.size(), ipoint / x_.size()};
}

size_t RegularGrid2D::cell_centered_grid_index_ip_jp(size_t i, size_t j) const {
    if (i >= x_.size() - 1 || j >= y_.size() - 1) {
        std::ostringstream oss;
        oss << "Invalid RegularGrid2d::cell_centered_grid_index_ip_jp() "
               "arguments: ";
        oss << std::endl << "    ";
        oss << "i=" << i << ", j=" << j;
        throw std::runtime_error(oss.str());
    }
    return i + j * (x_.size() - 1);
}

size_t RegularGrid2D::xface_grid_index_ip_j(size_t i, size_t j) const {
    if (i >= x_.size() - 1 || j >= y_.size()) {
        std::ostringstream oss;
        oss << "Invalid RegularGrid2d::xface_grid_index_ip_j() arguments: ";
        oss << std::endl << "    ";
        oss << "i=" << i << ", j=" << j;
        throw std::runtime_error(oss.str());
    }
    return i + j * (x_.size() - 1);
}

size_t RegularGrid2D::yface_grid_index_i_jp(size_t i, size_t j) const {
    if (i >= x_.size() || j >= y_.size() - 1) {
        std::ostringstream oss;
        oss << "Invalid RegularGrid2d::yface_grid_index_ip_j() arguments: ";
        oss << std::endl << "    ";
        oss << "i=" << i << ", j=" << j;
        throw std::runtime_error(oss.str());
    }
    return i + j * x_.size();
}

RegularGrid2D RegularGrid2D::cell_centered_grid() const {
    std::vector<double> ret_x;
    for (size_t i = 0; i < x_.size() - 1; ++i) {
        ret_x.push_back((x_[i] + x_[i + 1]) / 2);
    }
    std::vector<double> ret_y;
    for (size_t i = 0; i < y_.size() - 1; ++i) {
        ret_y.push_back((y_[i] + y_[i + 1]) / 2);
    }
    return RegularGrid2D(ret_x, ret_y);
}

RegularGrid2D RegularGrid2D::xface_centered_grid() const {
    std::vector<double> ret_x;
    for (size_t i = 0; i < x_.size() - 1; ++i) {
        ret_x.push_back((x_[i] + x_[i + 1]) / 2);
    }
    return RegularGrid2D(ret_x, y_);
}

RegularGrid2D RegularGrid2D::yface_centered_grid() const {
    std::vector<double> ret_y;
    for (size_t i = 0; i < y_.size() - 1; ++i) {
        ret_y.push_back((y_[i] + y_[i + 1]) / 2);
    }
    return RegularGrid2D(x_, ret_y);
}

bool RegularGrid2D::is_active_cell(size_t icell) const {
    return actnum_[icell];
}

void RegularGrid2D::deactivate_cells(Point bot_left, Point top_right) {
    std::vector<double> xcenters, ycenters;
    for (size_t i = 0; i < x_.size() - 1; ++i) {
        xcenters.push_back((x_[i] + x_[i + 1]) / 2);
    }
    for (size_t j = 0; j < y_.size() - 1; ++j) {
        ycenters.push_back((y_[j] + y_[j + 1]) / 2);
    }
    size_t i_begin = std::upper_bound(xcenters.begin(), xcenters.end(), bot_left.x) - xcenters.begin();
    size_t i_end = std::upper_bound(xcenters.begin(), xcenters.end(), top_right.x) - xcenters.begin();
    size_t j_begin = std::lower_bound(ycenters.begin(), ycenters.end(), bot_left.y) - ycenters.begin();
    size_t j_end = std::lower_bound(ycenters.begin(), ycenters.end(), top_right.y) - ycenters.begin();

    for (size_t i = i_begin; i < i_end; ++i)
        for (size_t j = j_begin; j < j_end; ++j) {
            actnum_[i + j * nx()] = 0;
        }
    set_face_types();
}

const std::vector<char>& RegularGrid2D::actnum() const {
    return actnum_;
}

RegularGrid2D::FaceType RegularGrid2D::yface_type(size_t yface_index) const {
    return yface_types_[yface_index];
}

RegularGrid2D::FaceType RegularGrid2D::xface_type(size_t xface_index) const {
    return xface_types_[xface_index];
}

const std::vector<RegularGrid2D::split_index_t>& RegularGrid2D::boundary_xfaces() const {
    return boundary_xfaces_;
}

const std::vector<RegularGrid2D::split_index_t>& RegularGrid2D::boundary_yfaces() const {
    return boundary_yfaces_;
}

void RegularGrid2D::set_face_types() {
    // ==== yfaces;
    yface_types_.resize((nx() + 1) * ny(), FaceType::Internal);
    boundary_yfaces_.clear();
    // left/right boundary
    for (size_t j = 0; j < ny(); ++j) {
        size_t ind_left = yface_grid_index_i_jp(0, j);
        size_t ind_right = yface_grid_index_i_jp(nx(), j);
        yface_types_[ind_left] = FaceType::Boundary;
        yface_types_[ind_right] = FaceType::Boundary;
        boundary_yfaces_.push_back({0, j});
        boundary_yfaces_.push_back({nx(), j});
    }
    // internal faces
    for (size_t i = 1; i < nx(); ++i)
        for (size_t j = 0; j < ny(); ++j) {
            size_t face_index = yface_grid_index_i_jp(i, j);
            size_t cell_left = cell_centered_grid_index_ip_jp(i - 1, j);
            size_t cell_right = cell_centered_grid_index_ip_jp(i, j);
            bool active_left = is_active_cell(cell_left);
            bool active_right = is_active_cell(cell_right);
            if (active_left && active_right) {
                // pass
            } else if (!active_left && !active_right) {
                yface_types_[face_index] = FaceType::Deactivated;
            } else {
                yface_types_[face_index] = FaceType::Boundary;
                boundary_yfaces_.push_back({i, j});
            }
        }
    // ==== xfaces;
    xface_types_.resize((ny() + 1) * nx(), FaceType::Internal);
    boundary_xfaces_.clear();
    // top/bot boundary
    for (size_t i = 0; i < nx(); ++i) {
        size_t ind_bot = xface_grid_index_ip_j(i, 0);
        size_t ind_top = xface_grid_index_ip_j(i, ny());
        xface_types_[ind_bot] = FaceType::Boundary;
        xface_types_[ind_top] = FaceType::Boundary;
        boundary_xfaces_.push_back({i, 0});
        boundary_xfaces_.push_back({i, ny()});
    }
    // internal faces
    for (size_t i = 0; i < nx(); ++i)
        for (size_t j = 1; j < ny(); ++j) {
            size_t face_index = xface_grid_index_ip_j(i, j);
            size_t cell_bot = cell_centered_grid_index_ip_jp(i, j - 1);
            size_t cell_top = cell_centered_grid_index_ip_jp(i, j);
            bool active_bot = is_active_cell(cell_bot);
            bool active_top = is_active_cell(cell_top);
            if (active_bot && active_top) {
                // pass
            } else if (!active_bot && !active_top) {
                xface_types_[face_index] = FaceType::Deactivated;
            } else {
                xface_types_[face_index] = FaceType::Boundary;
                boundary_xfaces_.push_back({i, j});
            }
        }
}
