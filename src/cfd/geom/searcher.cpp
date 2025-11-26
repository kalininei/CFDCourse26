#include "searcher.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <ranges>

using namespace cfd;
namespace bg = boost::geometry;

template<size_t Dim>
struct PointSearcher<Dim>::Impl {
    using point_t = bg::model::point<double, Dim, bg::cs::cartesian>;
    using value_t = std::pair<point_t, size_t>;

    static point_t build_point_t(Point p) {
        if constexpr (Dim == 1) {
            return point_t{p.x};
        } else if constexpr (Dim == 2) {
            return point_t{p.x, p.y};
        } else {
            return point_t{p.x, p.y, p.z};
        }
    }

    void add(const std::vector<Point>& points) {
        for (size_t i = 0; i < points.size(); ++i) {
            rtree_.insert(std::make_pair(build_point_t(points[i]), i));
        }
    }

    std::vector<size_t> nearest(const Point& p, size_t n) const {
        std::vector<value_t> vals;
        point_t point = build_point_t(p);
        rtree_.query(bg::index::nearest(point, n), std::back_inserter(vals));
        std::vector<size_t> ret;
        for (const auto& v: vals) {
            ret.push_back(v.second);
        }
        return ret;
    }

private:
    bg::index::rtree<value_t, bg::index::quadratic<16>> rtree_;
};

template<size_t Dim>
PointSearcher<Dim>::PointSearcher() {
    pimpl_ = std::make_shared<Impl>();
}

template<size_t Dim>
PointSearcher<Dim>::PointSearcher(const std::vector<Point>& points) {
    pimpl_ = std::make_shared<Impl>();
    add_points(points);
}

template<size_t Dim>
void PointSearcher<Dim>::add_points(const std::vector<Point>& points) {
    pimpl_->add(points);
}

template<size_t Dim>
std::vector<size_t> PointSearcher<Dim>::nearest(const Point& p, size_t n) const {
    return pimpl_->nearest(p, n);
}

template<size_t Dim>
struct BoxSearcher<Dim>::Impl {
    using point_t = bg::model::point<double, Dim, bg::cs::cartesian>;
    using box_t = bg::model::box<point_t>;
    using value_t = std::pair<box_t, size_t>;

    static point_t build_point_t(Point p, double margin = 0.0) {
        if constexpr (Dim == 1) {
            return point_t{p.x + margin};
        } else if constexpr (Dim == 2) {
            return point_t{p.x + margin, p.y + margin};
        } else {
            return point_t{p.x + margin, p.y + margin, p.z + margin};
        }
    }

    static box_t build_box_t(const std::array<Point, 2>& box) {
        return {build_point_t(box[0], -1e-12), build_point_t(box[1], 1e-12)};
    }

    void add(const std::vector<std::array<Point, 2>>& boxes) {
        for (size_t i = 0; i < boxes.size(); ++i) {
            rtree_.insert(std::make_pair(build_box_t(boxes[i]), i));
        }
    }

    std::vector<size_t> within(const Point& p) const {
        std::vector<value_t> vals;
        point_t point = build_point_t(p);
        rtree_.query(bg::index::contains(point), std::back_inserter(vals));
        auto ret = vals | std::views::values;
        return std::vector<size_t>(ret.begin(), ret.end());
    }

private:
    bg::index::rtree<value_t, bg::index::quadratic<16>> rtree_;
};

template<size_t Dim>
BoxSearcher<Dim>::BoxSearcher() {
    pimpl_ = std::make_shared<Impl>();
}

template<size_t Dim>
BoxSearcher<Dim>::BoxSearcher(const std::vector<std::array<Point, 2>>& boxes) {
    pimpl_ = std::make_shared<Impl>();
    add_boxes(boxes);
}

template<size_t Dim>
void BoxSearcher<Dim>::add_boxes(const std::vector<std::array<Point, 2>>& boxes) {
    pimpl_->add(boxes);
}

template<size_t Dim>
std::vector<size_t> BoxSearcher<Dim>::within(const Point& p) const {
    return pimpl_->within(p);
}

namespace cfd {
template struct PointSearcher<1>;
template struct PointSearcher<2>;
template struct PointSearcher<3>;
template struct BoxSearcher<1>;
template struct BoxSearcher<2>;
template struct BoxSearcher<3>;
} // namespace cfd
