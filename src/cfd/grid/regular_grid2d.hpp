#ifndef CFD_GEOM_REGULAR_GRID2_HPP
#define CFD_GEOM_REGULAR_GRID2_HPP

#include "cfd/grid/i_grid.hpp"

namespace cfd {

/**
 * @brief structured 2D quadraliteral grid with
 *
 * Fast indexing for points and cells goes along x dimension, so that the linear index equals
 * \f$
 * 	i = i_x + (n_x+1) i_y
 * \f$
 */
class RegularGrid2D : public IGrid2D {
public:
    /**
     * @brief Split index as the [ix, iy] pair
     */
    using split_index_t = std::array<size_t, 2>;

    /**
     * @brief Builds regular equidistant 2d grid
     *
     * @param x0  x minimum
     * @param x1  x maximum
     * @param y0  y minimum
     * @param y1  y maximum
     * @param nx  number of cells in x direction
     * @param ny  number of cells in y direction
     */
    RegularGrid2D(double x0, double x1, double y0, double y1, size_t nx, size_t ny);

    /**
     * @brief Builds non-equidistant regular 2d grid
     *
     * @param x  x-point values
     * @param y  y-point values
     *
     * The resulting number of points in the grid will be `x.size()*y.size()`
     */
    RegularGrid2D(const std::vector<double> &x, const std::vector<double> &y);

    double Lx() const;
    double Ly() const;

    /**
     * @brief Builds grid with its points located in the centers of this grid cells
     */
    RegularGrid2D cell_centered_grid() const;

    /**
     * @brief Builds grid with its points located in the centers of x-faces
     */
    RegularGrid2D xface_centered_grid() const;

    /**
     * @brief Builds grid with its points located in the centers of y-faces
     */
    RegularGrid2D yface_centered_grid() const;

    /**
     * @brief Converts split [ix, iy] point index to linear point index
     *
     * @param point_split_index  split point index as [ix, iy] array
     * @return linear point index
     */
    size_t to_linear_point_index(split_index_t point_split_index) const;

    /**
     * @brief Converts linear point index to [ix, iy] point index
     *
     * @param ipoint  linear point index in [0, n_points()-1]
     * @return split point index as the [ix, iy] array
     */
    split_index_t to_split_point_index(size_t ipoint) const;

    /**
     * @brief builds linear index for [i+1/2, j+1/2] point in the cell_centered grid.
     *
     * @param i  x-index
     * @param j  y-index
     *
     * Computes j*nx + i. Throws if i, j indices are not valid
     */
    size_t cell_centered_grid_index_ip_jp(size_t i, size_t j) const;

    /**
     * @brief builds linear index for [i+1/2, j] point in x-face grid.
     *
     * @param i  x-index
     * @param j  y-index
     *
     * Computes j*nx + i. Throws if i, j indices are not valid
     */
    size_t xface_grid_index_ip_j(size_t i, size_t j) const;

    /**
     * @brief builds linear index for [i, j+1/2] point in y-face grid.
     *
     * @param i  x-index
     * @param j  y-index
     *
     * Computes j*(nx+1) + i. Throws if i, j indices are not valid
     */
    size_t yface_grid_index_i_jp(size_t i, size_t j) const;

    /**
     * @brief Number of cells in x-direction
     */
    size_t nx() const;

    /**
     * @brief Number of cells in y-direction
     */
    size_t ny() const;

    /**
     * @brief returns acitve status of the cell
     */
    bool is_active_cell(size_t icell) const;

    /**
     * @brief deactivates all cells within bot_left, top_right rectangle
     */
    void deactivate_cells(Point bot_left, Point top_right);

    /**
     * @brief actnum vector for cells
     */
    const std::vector<char> &actnum() const;

    enum struct FaceType {
        Internal,
        Boundary,
        Deactivated
    };

    FaceType yface_type(size_t yface_index) const;
    FaceType xface_type(size_t xface_index) const;
    const std::vector<split_index_t> &boundary_yfaces() const;
    const std::vector<split_index_t> &boundary_xfaces() const;

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
    std::vector<double> _x;
    std::vector<double> _y;
    std::vector<char> _actnum;
    std::vector<FaceType> _yface_types;
    std::vector<FaceType> _xface_types;
    std::vector<split_index_t> _boundary_xfaces;
    std::vector<split_index_t> _boundary_yfaces;

    void set_face_types();
};

} // namespace cfd
#endif
