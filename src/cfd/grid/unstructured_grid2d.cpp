#include "unstructured_grid2d.hpp"
#include "cfd/geom/simplex.hpp"
#include "cfd/grid/vtk.hpp"
#include <cassert>
#include <fstream>
#include <sstream>

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Cache
///////////////////////////////////////////////////////////////////////////////

void UnstructuredGrid2D::Cache::clear() {
    cell_centers.clear();
    cell_volumes.clear();
    face_normals.clear();
    face_areas.clear();
    tab_cell_face.clear();
}

void UnstructuredGrid2D::Cache::need_cell_centers(const UnstructuredGrid2D& grid) {
    if (cell_centers.size() > 0) {
        return;
    }
    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        double sum_area = 0;
        double sum_x = 0;
        double sum_y = 0;
        const auto& cp = grid.tab_cell_point(icell);
        Point p0 = grid.point(cp[0]);
        for (size_t i = 1; i < cp.size() - 1; ++i) {
            Point p1 = grid.point(cp[i]);
            Point p2 = grid.point(cp[i + 1]);
            double area = triangle_area(p0, p1, p2);
            double x = (p0.x() + p1.x() + p2.x()) / 3.0;
            double y = (p0.y() + p1.y() + p2.y()) / 3.0;
            sum_x += x * area;
            sum_y += y * area;
            sum_area += area;
        }
        cell_centers.push_back({sum_x / sum_area, sum_y / sum_area});
        cell_volumes.push_back(sum_area);
    }
}

void UnstructuredGrid2D::Cache::need_cell_volumes(const UnstructuredGrid2D& grid) {
    return need_cell_centers(grid);
}

void UnstructuredGrid2D::Cache::need_face_normals(const UnstructuredGrid2D& grid) {
    if (face_normals.size() > 0) {
        return;
    }
    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        Point p0 = grid.point(grid._face_points[iface][0]);
        Point p1 = grid.point(grid._face_points[iface][1]);
        Vector s = p1 - p0;
        face_normals.push_back(Vector(s.y(), -s.x()) / vector_abs(s));
    }
}

void UnstructuredGrid2D::Cache::need_face_areas(const UnstructuredGrid2D& grid) {
    if (face_areas.size() > 0) {
        return;
    }
    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        Point p0 = grid.point(grid._face_points[iface][0]);
        Point p1 = grid.point(grid._face_points[iface][1]);
        double d = vector_abs(p1 - p0);
        face_areas.push_back(d);
    }
}

void UnstructuredGrid2D::Cache::need_tab_cell_face(const UnstructuredGrid2D& grid) {
    if (tab_cell_face.size() > 0) {
        return;
    }
    tab_cell_face.resize(grid.n_cells());
    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        std::array<size_t, 2> cells = grid.tab_face_cell(iface);
        if (cells[0] != INVALID_INDEX) {
            tab_cell_face[cells[0]].push_back(iface);
        }
        if (cells[1] != INVALID_INDEX) {
            tab_cell_face[cells[1]].push_back(iface);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Unstructured Grid
///////////////////////////////////////////////////////////////////////////////

UnstructuredGrid2D::UnstructuredGrid2D(const std::vector<Point>& points,
                                       const std::vector<std::vector<size_t>>& cell_point)
    : _points(points), _cells(cell_point) {
    initialize();
}

namespace {

std::vector<std::vector<size_t>> assemble_cell_point(const IGrid2D& grid) {
    std::vector<std::vector<size_t>> ret;
    for (size_t i = 0; i < grid.n_cells(); ++i) {
        ret.push_back(grid.tab_cell_point(i));
    }
    return ret;
}

} // namespace

UnstructuredGrid2D::UnstructuredGrid2D(const IGrid2D& grid)
    : UnstructuredGrid2D(grid.points(), assemble_cell_point(grid)) {}

void UnstructuredGrid2D::validate() {
    if (_cells.size() == 0) {
        throw std::runtime_error("no cells in grid");
    }
    if (_points.size() == 0) {
        throw std::runtime_error("no points in grid");
    }
}

void UnstructuredGrid2D::initialize() {
    _cache.clear();
    validate();
    // assemble faces
    std::vector<std::array<size_t, 3>> point_point_cell;
    for (size_t icell = 0; icell < n_cells(); ++icell) {
        size_t p0 = _cells[icell].back();
        for (size_t i = 0; i < _cells[icell].size(); ++i) {
            size_t p1 = _cells[icell][i];
            point_point_cell.push_back({p0, p1, icell});
            std::swap(p0, p1);
        }
    }
    auto sorter_less = [](const std::array<size_t, 3>& v1, const std::array<size_t, 3>& v2) -> bool {
        size_t min1 = std::min(v1[0], v1[1]);
        size_t max1 = std::max(v1[0], v1[1]);
        size_t min2 = std::min(v2[0], v2[1]);
        size_t max2 = std::max(v2[0], v2[1]);
        if (min1 < min2) {
            return true;
        } else if (min1 > min2) {
            return false;
        } else if (max1 < max2) {
            return true;
        } else {
            return false;
        }
    };
    std::sort(point_point_cell.begin(), point_point_cell.end(), sorter_less);

    auto add_face = [this](const std::array<size_t, 3>& v) {
        _face_points.push_back({v[0], v[1]});
        _face_cells.push_back({v[2], INVALID_INDEX});
    };
    for (size_t i = 0; i < point_point_cell.size(); ++i) {
        if (i == 0 || sorter_less(point_point_cell[i - 1], point_point_cell[i])) {
            add_face(point_point_cell[i]);
        } else {
            _face_cells.back()[1] = point_point_cell[i][2];
        }
    }
    for (size_t i = 0; i < _face_points.size(); ++i) {
        if (_face_points[i][0] > _face_points[i][1]) {
            std::swap(_face_points[i][0], _face_points[i][1]);
            std::swap(_face_cells[i][0], _face_cells[i][1]);
        }
    }
}

size_t UnstructuredGrid2D::n_points() const {
    return _points.size();
}

size_t UnstructuredGrid2D::n_cells() const {
    return _cells.size();
}

size_t UnstructuredGrid2D::n_faces() const {
    return _face_cells.size();
}

Point UnstructuredGrid2D::point(size_t ipoint) const {
    return _points[ipoint];
}

Point UnstructuredGrid2D::cell_center(size_t icell) const {
    _cache.need_cell_centers(*this);
    return _cache.cell_centers[icell];
}

double UnstructuredGrid2D::cell_volume(size_t icell) const {
    _cache.need_cell_volumes(*this);
    return _cache.cell_volumes[icell];
}

Vector UnstructuredGrid2D::face_normal(size_t iface) const {
    _cache.need_face_normals(*this);
    return _cache.face_normals[iface];
}

double UnstructuredGrid2D::face_area(size_t iface) const {
    _cache.need_face_areas(*this);
    return _cache.face_areas[iface];
}

Point UnstructuredGrid2D::face_center(size_t iface) const {
    Point p0 = point(_face_points[iface][0]);
    Point p1 = point(_face_points[iface][1]);
    return (p0 + p1) / 2.0;
}

std::vector<Point> UnstructuredGrid2D::points() const {
    return _points;
}

std::vector<size_t> UnstructuredGrid2D::tab_cell_point(size_t icell) const {
    return _cells[icell];
}

std::array<size_t, 2> UnstructuredGrid2D::tab_face_cell(size_t iface) const {
    return _face_cells[iface];
}

std::vector<size_t> UnstructuredGrid2D::tab_face_point(size_t iface) const {
    return std::vector<size_t>(_face_points[iface].begin(), _face_points[iface].end());
}

std::vector<size_t> UnstructuredGrid2D::tab_cell_face(size_t icell) const {
    _cache.need_tab_cell_face(*this);
    return _cache.tab_cell_face[icell];
}

void UnstructuredGrid2D::save_vtk(std::string fname) const {
    std::ofstream fs(fname);
    VtkUtils::append_header("Grid2", fs);
    VtkUtils::append_points(this->points(), fs);
    // Cells
    std::ostringstream oss;
    size_t n_totals = 0;
    for (size_t i = 0; i < this->n_cells(); ++i) {
        auto v = this->tab_cell_point(i);
        oss << v.size();
        for (size_t p : v) {
            oss << " " << p;
        }
        oss << std::endl;
        n_totals += v.size() + 1;
    }
    fs << "CELLS " << this->n_cells() << " " << n_totals << std::endl;
    fs << oss.str();
    fs << "CELL_TYPES " << this->n_cells() << std::endl;
    for (size_t i = 0; i < this->n_cells(); ++i)
        fs << 7 << std::endl;
}

namespace {

struct ELineNotFound : public std::runtime_error {
    ELineNotFound(std::string s) : std::runtime_error(s + " line not found while reading input") {};
};

std::string get_line_by_start(std::string start, std::istream& is) {
    std::string line;
    while (is) {
        std::getline(is, line);
        if (!line.empty() && line.back() == '\r') {
            line.resize(line.size() - 1);
        }
        if (line.substr(0, start.size()) == start) {
            return line;
        }
    }
    throw ELineNotFound(start);
}

} // namespace
UnstructuredGrid2D UnstructuredGrid2D::vtk_read(std::string filename, bool silent) {
    if (!silent)
        std::cout << "Reading grid from " << filename << std::endl;

    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error(filename + " is not found");
    }
    std::string line, tmp;
    // header
    line = get_line_by_start("DATASET", ifs);
    if (line.substr(8) != "UNSTRUCTURED_GRID") {
        throw std::runtime_error("Only unstructured grid can be read");
    }
    // points
    line = get_line_by_start("POINTS", ifs);
    int n_points;
    std::istringstream(line) >> tmp >> n_points;

    std::vector<Point> points(n_points);
    double z;
    for (int i = 0; i < n_points; ++i) {
        ifs >> points[i].x() >> points[i].y() >> z;
        if (std::abs(z) > 0) {
            throw std::runtime_error("Z-coordinate for 2d grids should be zero");
        }
    }
    if (!silent)
        std::cout << "-- " << n_points << " points" << std::endl;

    // cells
    int n_totals, n_cells;
    line = get_line_by_start("CELLS", ifs);
    std::istringstream(line) >> tmp >> n_cells >> n_totals;
    std::vector<int> totals(n_totals);
    std::vector<int> types(n_cells);
    // -- totals
    for (int i = 0; i < n_totals; ++i) {
        ifs >> totals[i];
    }
    // -- types
    line = get_line_by_start("CELL_TYPES", ifs);
    for (int i = 0; i < n_cells; ++i) {
        ifs >> types[i];
    }

    std::vector<std::vector<size_t>> cell_points;
    int* cursor = &totals[0];
    // assemble vtk cells
    for (int i = 0; i < n_cells; ++i) {
        int len = cursor[0];
        cursor++;
        int type = types[i];
        // std::vector<int> points(cursor, cursor + len);

        switch (type) {
        case 5:
        case 7:
        case 9:
            cell_points.emplace_back(cursor, cursor + len);
            break;
        case 6:
            throw std::runtime_error("Triangle strips are not supported");
        case 8:
            cell_points.push_back(
                {(size_t)*cursor, (size_t)*(cursor + 1), (size_t)*(cursor + 3), (size_t)*(cursor + 2)});
            break;
        default:
            break;
        }

        cursor += len;
    }
    if (!silent)
        std::cout << "-- " << cell_points.size() << " 2d cells" << std::endl;

    UnstructuredGrid2D ret(points, cell_points);

    // check cell direction
    int bad_dir = 0;
    for (size_t i = 0; i < ret.n_cells(); ++i) {
        double v = ret.cell_volume(i);
        if (v < 0) {
            bad_dir++;
            std::reverse(cell_points[i].begin() + 1, cell_points[i].end());
        }
    }
    if (bad_dir > 0) {
        if (!silent)
            std::cout << "-- " << bad_dir << " cells are reversed" << std::endl;
        return UnstructuredGrid2D(points, cell_points);
    } else {
        return ret;
    }
}
