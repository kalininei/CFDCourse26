#include "polygon.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>

using namespace cfd;

namespace boost::geometry::traits {
template<>
struct point_order<std::vector<Point>> {
    static const order_selector value = counterclockwise;
};

template<>
struct closure<std::vector<Point>> {
    static const closure_selector value = open;
};
} // namespace boost::geometry::traits

BOOST_GEOMETRY_REGISTER_POINT_2D(Point, double, boost::geometry::cs::cartesian, x, y)
BOOST_GEOMETRY_REGISTER_RING(std::vector<Point>)

bool cfd::is_in_polygon(Point p, const std::vector<Point>& polygon, double eps) {
    return polygon.size() > 2 && (boost::geometry::within(p, polygon) || boost::geometry::distance(p, polygon) <= eps);
}
