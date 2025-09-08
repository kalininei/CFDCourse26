#ifndef CFD_GRID_I_GRID_HPP
#define CFD_GRID_I_GRID_HPP

#include "cfd/cfd_common.hpp"
#include "cfd/geom/point.hpp"

namespace cfd {

/**
 * @brief Any dimension grid interface
 */
class IGrid {
public:
    virtual ~IGrid() = default;

    /// @brief number of geometric dimensions: 1, 2 or 3
    virtual size_t dim() const = 0;

    /// @brief number of grid points
    virtual size_t n_points() const = 0;

    /// @brief number of grid cells
    virtual size_t n_cells() const = 0;

    /// @brief number of grid faces
    virtual size_t n_faces() const = 0;

    /// @brief get point at i-th index
    virtual Point point(size_t ipoint) const = 0;

    /// @brief get cell center
    virtual Point cell_center(size_t icell) const = 0;

    /// @brief cell volume
    virtual double cell_volume(size_t icell) const = 0;

    /// @brief linear cell size
    virtual double cell_size(size_t icell) const = 0;

    /// @brief get face unit normal
    virtual Vector face_normal(size_t iface) const = 0;

    /// @brief get face area
    virtual double face_area(size_t iface) const = 0;

    /// @brief get face center
    virtual Point face_center(size_t iface) const = 0;

    /// @brief get list of all points
    virtual std::vector<Point> points() const = 0;

    /**
     * @brief cell->point connectivity table
     *
     * @param icell  cell index
     *
     * For 1d grid points are ordered by x coordinate \n
     * For 2d grid points have counter clockwise direction\n
     * 3d grids have no special points ordering
     */
    virtual std::vector<size_t> tab_cell_point(size_t icell) const = 0;

    /**
     * @brief face->cell connectivity table
     *
     * @param iface  face index
     *
     * @return pair of cell indexes adjacent to given cell
     *
     * The first returned cell index corresponds to negative face side, and
     * the second - to positive. Positive face side lies towards face normal.
     *
     * If this is a boundary face then the respective cell index will be equal to
     * INVALID_INDEX.
     */
    virtual std::array<size_t, 2> tab_face_cell(size_t iface) const = 0;

    virtual std::vector<size_t> tab_face_point(size_t iface) const = 0;
    virtual std::vector<size_t> tab_cell_face(size_t icell) const = 0;

    virtual std::vector<size_t> boundary_faces() const;
    virtual std::vector<size_t> boundary_points() const;

    /**
     * @brief Saves grid to vtk format
     *
     * @param fname  output filename
     */
    virtual void save_vtk(std::string fname) const = 0;

private:
    struct Cache {
        std::vector<size_t> boundary_faces;
        std::vector<size_t> boundary_points;

        void clear();
        void need_boundary_faces(const IGrid& grid);
        void need_boundary_points(const IGrid& grid);
    };
    mutable Cache _cache;
};

/**
 * @brief Abstract 1D grid
 */
class IGrid1D : public IGrid {
public:
    virtual ~IGrid1D() = default;

    size_t dim() const override {
        return 1;
    }

    double cell_size(size_t icell) const override {
        return cell_volume(icell);
    }
};

/**
 * @brief Abstract 2D grid
 */
class IGrid2D : public IGrid {
public:
    virtual ~IGrid2D() = default;

    size_t dim() const override {
        return 2;
    }

    double cell_size(size_t icell) const override {
        return std::sqrt(cell_volume(icell));
    }
};

/**
 * @brief Abstract 3D grid
 */
class IGrid3D : public IGrid {
public:
    virtual ~IGrid3D() = default;

    size_t dim() const override {
        return 3;
    }

    double cell_size(size_t icell) const override {
        return std::cbrt(cell_volume(icell));
    }
};

} // namespace cfd
#endif
