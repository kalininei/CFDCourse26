#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/regular_grid2d.hpp"
#include "cfd/grid/unstructured_grid2d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd26_test.hpp"
#include "utils/filesystem.hpp"

using namespace cfd;

TEST_CASE("Grid1d", "[grid1]") {
    Grid1D grid(0, 1, 3);

    CHECK(grid.n_points() == 4);
    CHECK(grid.n_cells() == 3);
    CHECK(grid.n_faces() == 4);

    grid.save_vtk("out1.vtk");
    VtkUtils::add_point_data(std::vector<double>{0, 1, 2, 3}, "data1", "out1.vtk");
}

TEST_CASE("RegularGrid2d", "[reggrid2]") {
    RegularGrid2D grid(0, 1, 1, 3, 3, 2);

    CHECK(grid.n_points() == 12);
    CHECK(grid.n_cells() == 6);

    CHECK(grid.to_split_point_index(8)[0] == 0);
    CHECK(grid.to_split_point_index(8)[1] == 2);
    CHECK(grid.to_split_point_index(11)[0] == 3);
    CHECK(grid.to_split_point_index(11)[1] == 2);
    CHECK(grid.to_linear_point_index({0, 0}) == 0);
    CHECK(grid.to_linear_point_index({2, 0}) == 2);
    CHECK(grid.to_linear_point_index({3, 1}) == 7);
    CHECK(grid.cell_center(0).x() == Approx(0.1666).margin(1e-2));
    CHECK(grid.cell_center(0).y() == Approx(1.5).margin(1e-2));
    CHECK(grid.cell_center(4).x() == Approx(0.5).margin(1e-2));
    CHECK(grid.cell_center(4).y() == Approx(2.5).margin(1e-2));

    grid.save_vtk("out2.vtk");
    VtkUtils::add_point_data(std::vector<double>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, "data1", "out2.vtk");
}

TEST_CASE("RegularGrid2d-nonuniform", "[reggrid2-nonuni]") {
    std::vector<double> x{10, 12, 15, 16};
    std::vector<double> y{-10, -8, -3};
    RegularGrid2D grid(x, y);

    CHECK(grid.Lx() == Approx(6));
    CHECK(grid.Ly() == Approx(7));
    CHECK(grid.cell_center(4).x() == Approx(13.5));
    CHECK(grid.cell_center(4).y() == Approx(-5.5));

    CHECK(grid.n_faces() == 17);
    CHECK(grid.face_center(8).x() == Approx(15.5));
    CHECK(grid.face_center(8).y() == Approx(-3));
    CHECK(grid.face_center(13).x() == Approx(10));
    CHECK(grid.face_center(13).y() == Approx(-5.5));
    CHECK(grid.face_center(16).x() == Approx(16));
    CHECK(grid.face_center(16).y() == Approx(-5.5));
    CHECK(grid.face_normal(4).x() == Approx(0.0));
    CHECK(grid.face_normal(4).y() == Approx(1.0));
    CHECK(grid.face_normal(10).x() == Approx(1.0));
    CHECK(grid.face_normal(10).y() == Approx(0.0));
    CHECK(grid.face_area(6) == Approx(2));
    CHECK(grid.face_area(12) == Approx(2));
}

TEST_CASE("UnstructuredGrid2d", "[unstructured-grid2]") {
    std::vector<double> x{10, 12, 15, 16};
    std::vector<double> y{-10, -8, -3};
    RegularGrid2D grid(x, y);
    UnstructuredGrid2D ugrid(grid);

    CHECK(ugrid.n_cells() == 6);
    CHECK(ugrid.n_points() == 12);
    CHECK(ugrid.n_faces() == 17);

    CHECK(ugrid.cell_center(4).x() == Approx(13.5));
    CHECK(ugrid.cell_center(4).y() == Approx(-5.5));

    CHECK(ugrid.face_center(16).x() == Approx(15.5));
    CHECK(ugrid.face_center(16).y() == Approx(-3));
    CHECK(ugrid.face_center(8).x() == Approx(10));
    CHECK(ugrid.face_center(8).y() == Approx(-5.5));
    CHECK(ugrid.face_center(13).x() == Approx(16));
    CHECK(ugrid.face_center(13).y() == Approx(-5.5));
    CHECK(ugrid.face_normal(9).x() == Approx(0.0));
    CHECK(ugrid.face_normal(9).y() == Approx(-1.0));
    CHECK(ugrid.face_normal(3).x() == Approx(1.0));
    CHECK(ugrid.face_normal(3).y() == Approx(0.0));
    CHECK(ugrid.face_area(14) == Approx(2));
    CHECK(ugrid.face_area(6) == Approx(2));
}

TEST_CASE("UnstructuredGrid2d, read from vtk", "[unstructured2-vtk]") {
    std::string fn = test_directory_file("hexagrid_50.vtk");
    UnstructuredGrid2D grid = UnstructuredGrid2D::vtk_read(fn);

    CHECK(grid.n_points() == 106);
    CHECK(grid.n_cells() == 52);
    CHECK(grid.cell_volume(0) == Approx(0.00595238).margin(1e-6));
    CHECK(grid.cell_volume(34) == Approx(0.0238095).margin(1e-6));
    CHECK(grid.tab_face_cell(53)[0] == 5);
    CHECK(grid.tab_face_cell(53)[1] == 21);
    CHECK(grid.cell_center(36).x() == Approx(0.714286).margin(1e-6));
    CHECK(grid.cell_center(36).y() == Approx(0.583333).margin(1e-6));
    CHECK(grid.cell_center(0).x() == Approx(0.037037037).margin(1e-6));
    CHECK(grid.cell_center(0).y() == Approx(0.037037037).margin(1e-6));
    CHECK(grid.cell_center(7).y() == Approx((grid.point(12).y() + grid.point(13).y()) / 2.0).margin(1e-6));

    CHECK_NOTHROW(grid.save_vtk("hexa.vtk"));
}

TEST_CASE("Load grid from windows build", "[unstructured2-win]") {
    auto grid = UnstructuredGrid2D::vtk_read(test_directory_file("pebigrid_from_win.vtk"));
    CHECK(grid.n_cells() == 1022);
}
