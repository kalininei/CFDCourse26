#ifndef CFD_GEOM_SEARCHER_HPP
#define CFD_GEOM_SEARCHER_HPP

#include <memory>
#include <vector>
#include "cfd/geom/point.hpp"

namespace cfd{

template<size_t Dim=3>
class PointSearcher{
public:
	PointSearcher();
	PointSearcher(const std::vector<Point>& points);
	void add_points(const std::vector<Point>& points);

	std::vector<size_t> nearest(const Point& p, size_t n) const;
private:
	struct Impl;
	std::shared_ptr<Impl> _pimpl;
};

}
#endif
