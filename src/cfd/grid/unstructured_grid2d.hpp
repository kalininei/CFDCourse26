#ifndef CFD_UNSTRUCTURED_GRID_2D_HPP
#define CFD_UNSTRUCTURED_GRID_2D_HPP


#include "cfd/grid/i_grid.hpp"

namespace cfd{

class UnstructuredGrid2D: public IGrid2D{
public:
	/**
	 * @brief Builds unstructured 2d grid using points table and cell_point connectivity
	 *
	 * @param points      points vector
	 * @param cell_point  cell_point connectivity in counter_clockwise direction 
	 */
	UnstructuredGrid2D(const std::vector<Point>& points, const std::vector<std::vector<size_t>>& cell_point);

	/**
	 * @brief converts abstract 2d grid into unstructured format
	 */
	UnstructuredGrid2D(const IGrid2D& grid);

	/*
	 * @brief Read grid from vtk file
	 *
	 * @param filename  vtk ascii legacy format file name
	 */
	static UnstructuredGrid2D vtk_read(std::string filename, bool silent=false);

	// overridden methods
	size_t n_points() const override;
	size_t n_cells() const override;
	size_t n_faces() const override;
	Point point(size_t ipoint) const override;
	Point cell_center(size_t icell) const override;
	double cell_volume(size_t icell) const override;
	Vector face_normal(size_t iface) const override;
	double face_area(size_t iface) const override;
	Point face_center(size_t iface) const override;
	std::vector<Point> points() const override;
	std::vector<size_t> tab_cell_point(size_t icell) const override;
	std::array<size_t, 2> tab_face_cell(size_t iface) const override;
	std::vector<size_t> tab_face_point(size_t iface) const override;
	std::vector<size_t> tab_cell_face(size_t icell) const override;
	void save_vtk(std::string fname) const override;
private:
	const std::vector<Point> _points;
	const std::vector<std::vector<size_t>> _cells;
	std::vector<std::array<size_t, 2>> _face_cells;
	std::vector<std::array<size_t, 2>> _face_points;

	struct Cache{
		std::vector<Point> cell_centers;
		std::vector<double> cell_volumes;
		std::vector<Vector> face_normals;
		std::vector<double> face_areas;
		std::vector<std::vector<size_t>> tab_cell_face;

		void clear();
		void need_cell_centers(const UnstructuredGrid2D& grid);
		void need_cell_volumes(const UnstructuredGrid2D& grid);
		void need_face_normals(const UnstructuredGrid2D& grid);
		void need_face_areas(const UnstructuredGrid2D& grid);
		void need_tab_cell_face(const UnstructuredGrid2D& grid);
	};
	mutable Cache _cache;

	void validate();
	void initialize();
};

}
#endif
