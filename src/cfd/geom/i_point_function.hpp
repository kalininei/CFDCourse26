#ifndef __CFD_GEOM_POINT_FUNCTION_HPP__
#define __CFD_GEOM_POINT_FUNCTION_HPP__

#include "cfd/geom/point.hpp"

namespace cfd {

struct IPointFunction {
    virtual ~IPointFunction() = default;

    virtual double value(Point) const = 0;
    virtual Vector grad(Point) const = 0;

    double operator()(Point p) const {
        return value(p);
    };
};

} // namespace cfd

#endif
