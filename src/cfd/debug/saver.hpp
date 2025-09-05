#ifndef __CFD_DEBUG_SAVER_HPP__
#define __CFD_DEBUG_SAVER_HPP__

#include "cfd/grid/i_grid.hpp"
#include <vector>

namespace cfd {
namespace dbg {

void ping_saver_cpp();
void save_point_data(const IGrid &grid, const std::vector<double> &data);
void save_cell_data(const IGrid &grid, const std::vector<double> &data);
void save_cell_vector(const IGrid &grid, const std::vector<Vector> &vec);
void save_extended_colloc_data(const IGrid &grid, const std::vector<double> &data);
void save_extended_colloc_data(const IGrid &grid, const double *data);
void save_face_data(const IGrid &grid, const std::vector<double> &data);

} // namespace dbg
} // namespace cfd
#endif
