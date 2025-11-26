#ifndef __CFD_GRID_CELL_FINDER_HPP__
#define __CFD_GRID_CELL_FINDER_HPP__

#include "cfd/grid/i_grid.hpp"

namespace cfd {

class CellFinder {
public:
    CellFinder(const IGrid& grid);
    ~CellFinder();
    size_t operator()(Point p) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}; // namespace cfd

#endif
