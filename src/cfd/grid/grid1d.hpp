#ifndef CFD_GEOM_GRID1_HPP
#define CFD_GEOM_GRID1_HPP

#include "cfd/grid/i_grid.hpp"

namespace cfd{

/**
 * @brief 1D grid with ordered points
 */
class Grid1D: public IGrid1D{
public:
	/**
	 * Builds regular equidistant 1d grid
	 *
	 * @param left     left boundary
	 * @param right    right boundary
	 * @param n_cells  number of cells
	 */
	Grid1D(double left, double right, size_t n_cells);

	// overridden methods
	size_t n_points() const override;
	size_t n_cells() const override;
	size_t n_faces() const override;
	Point cell_center(size_t icell) const override;
	double cell_volume(size_t icell) const override;
	Vector face_normal(size_t iface) const override;
	double face_area(size_t iface) const override;
	Point face_center(size_t iface) const override;
	Point point(size_t ipoint) const override;
	std::vector<Point> points() const override;
	std::vector<size_t> tab_cell_point(size_t icell) const override;
	std::array<size_t, 2> tab_face_cell(size_t iface) const override;
	std::vector<size_t> tab_face_point(size_t iface) const override;
	std::vector<size_t> tab_cell_face(size_t icell) const override;
	void save_vtk(std::string fname) const override;
private:
	std::vector<Point> _points;
};

}
#endif
