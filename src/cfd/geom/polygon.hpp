#ifndef __CFD_GEOM_POLYGON_HPP__
#define __CFD_GEOM_POLYGON_HPP__

#include "cfd/geom/point.hpp"
#include <vector>

namespace cfd {

bool is_in_polygon(Point p, const std::vector<Point>& polygon, double eps = 1e-12);

};

#endif
