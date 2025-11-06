#include "saver.hpp"
#include "cfd/debug/tictoc.hpp"
#include "cfd/grid/vtk.hpp"

using namespace cfd;

namespace {

const std::string dbg_vtk_filename = "dbg.vtk";
const std::string dbg_vtk_cap = "data";

void save_faces(const IGrid& grid) {
    std::vector<Point> points;
    for (size_t i = 0; i < grid.n_points(); ++i) {
        points.push_back(grid.point(i));
    }

    std::ofstream fs(dbg_vtk_filename);
    VtkUtils::append_header("Grid2", fs);
    VtkUtils::append_points(points, fs);

    // Faces
    std::vector<int> types;
    std::ostringstream oss;
    size_t n_totals = 0;
    for (size_t i = 0; i < grid.n_faces(); ++i) {
        std::vector<size_t> pp = grid.tab_face_point(i);
        oss << pp.size();
        for (size_t p: pp) {
            oss << " " << p;
        }
        if (pp.size() == 1) {
            types.push_back(1);
        } else if (pp.size() == 2) {
            types.push_back(3);
        } else {
            types.push_back(7);
        }
        oss << std::endl;
        n_totals += pp.size() + 1;
    }
    fs << "CELLS " << grid.n_faces() << " " << n_totals << std::endl;
    fs << oss.str();
    fs << "CELL_TYPES " << grid.n_faces() << std::endl;

    for (size_t i = 0; i < grid.n_faces(); ++i)
        fs << types[i] << std::endl;
}

size_t save_extended_colloc(const IGrid& grid) {
    std::vector<Point> points;
    for (size_t i = 0; i < grid.n_points(); ++i) {
        points.push_back(grid.point(i));
    }

    std::ofstream fs(dbg_vtk_filename);
    VtkUtils::append_header("Grid2", fs);
    VtkUtils::append_points(points, fs);

    // Cells
    std::vector<int> types;
    std::ostringstream oss;
    size_t n_totals = 0;
    for (size_t i = 0; i < grid.n_cells(); ++i) {
        auto v = grid.tab_cell_point(i);
        oss << v.size();
        for (size_t p: v) {
            oss << " " << p;
        }
        oss << std::endl;
        n_totals += v.size() + 1;
        types.push_back(7);
    }
    // Faces
    for (size_t i = 0; i < grid.n_faces(); ++i) {
        std::array<size_t, 2> c = grid.tab_face_cell(i);
        if (c[0] != INVALID_INDEX && c[1] != INVALID_INDEX) {
            continue;
        }
        std::vector<size_t> pp = grid.tab_face_point(i);
        oss << pp.size();
        for (size_t p: pp) {
            oss << " " << p;
        }
        if (pp.size() == 1) {
            types.push_back(1);
        } else if (pp.size() == 2) {
            types.push_back(3);
        } else {
            types.push_back(7);
        }
        oss << std::endl;
        n_totals += pp.size() + 1;
    }
    fs << "CELLS " << types.size() << " " << n_totals << std::endl;
    fs << oss.str();
    fs << "CELL_TYPES " << types.size() << std::endl;

    for (size_t i = 0; i < types.size(); ++i) {
        fs << types[i] << std::endl;
    }

    return types.size();
}

} // namespace

void cfd::dbg::save_point_data(const IGrid& grid, const std::vector<double>& data) {
    grid.save_vtk(dbg_vtk_filename);
    std::vector<double> d2(data.begin(), data.begin() + grid.n_points());
    VtkUtils::add_point_data(d2, dbg_vtk_cap, dbg_vtk_filename);
}

void cfd::dbg::save_cell_data(const IGrid& grid, const std::vector<double>& data) {
    grid.save_vtk(dbg_vtk_filename);
    std::vector<double> d2(data.begin(), data.begin() + grid.n_cells());
    VtkUtils::add_cell_data(d2, dbg_vtk_cap, dbg_vtk_filename);
}

void cfd::dbg::save_extended_colloc_data(const IGrid& grid, const std::vector<double>& data) {
    size_t n_colloc = save_extended_colloc(grid);
    std::vector<double> d2(data.begin(), data.begin() + n_colloc);
    VtkUtils::add_cell_data(d2, dbg_vtk_cap, dbg_vtk_filename);
}

void cfd::dbg::save_extended_colloc_data(const IGrid& grid, const double* data) {
    size_t n_colloc = save_extended_colloc(grid);
    std::vector<double> d2(data, data + n_colloc);
    VtkUtils::add_cell_data(d2, dbg_vtk_cap, dbg_vtk_filename);
}

void cfd::dbg::save_face_data(const IGrid& grid, const std::vector<double>& data) {
    save_faces(grid);
    std::vector<double> d2(data.begin(), data.begin() + grid.n_faces());
    VtkUtils::add_cell_data(d2, dbg_vtk_cap, dbg_vtk_filename);
}

void cfd::dbg::save_cell_vector(const IGrid& grid, const std::vector<Vector>& vec) {
    grid.save_vtk(dbg_vtk_filename);
    std::vector<Vector> d2(vec.begin(), vec.begin() + grid.n_cells());
    VtkUtils::add_cell_vector(d2, dbg_vtk_cap, dbg_vtk_filename);
}

void cfd::dbg::ping_saver_cpp() {}
