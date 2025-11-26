#ifndef CFD_GEOM_SEARCHER_HPP
#define CFD_GEOM_SEARCHER_HPP

#include "cfd/geom/point.hpp"
#include <memory>
#include <vector>

namespace cfd {

template<size_t Dim = 3>
class PointSearcher {
    static_assert(Dim == 1 || Dim == 2 || Dim == 3);

public:
    PointSearcher();
    PointSearcher(const std::vector<Point>& points);
    void add_points(const std::vector<Point>& points);

    std::vector<size_t> nearest(const Point& p, size_t n) const;

private:
    struct Impl;
    std::shared_ptr<Impl> pimpl_;
};

template<size_t Dim = 3>
class BoxSearcher {
    static_assert(Dim == 1 || Dim == 2 || Dim == 3);

public:
    BoxSearcher();
    BoxSearcher(const std::vector<std::array<Point, 2>>& boxes);
    void add_boxes(const std::vector<std::array<Point, 2>>& boxes);

    std::vector<size_t> within(const Point& p) const;

private:
    struct Impl;
    std::shared_ptr<Impl> pimpl_;
};

} // namespace cfd
#endif
