#include "cell_finder.hpp"
#include "cfd/geom/polygon.hpp"
#include "cfd/geom/searcher.hpp"
#include <ranges>

using namespace cfd;

namespace {
constexpr double GEOM_EPS = 1e-12;
};

struct CellFinder::Impl {
    virtual ~Impl() = default;
    virtual size_t find_cell(Point p) const = 0;

    struct Finder1;
    struct Finder2;
};

/////////////////////////////////////////////////////////////////////
// 1D
/////////////////////////////////////////////////////////////////////
struct CellFinder::Impl::Finder1 : CellFinder::Impl {
    Finder1(const IGrid1D& g) {
        for (size_t i = 0; i < g.n_points(); ++i) {
            xvals_.push_back(g.point(i).x);
        }
    }

    size_t find_cell(Point p) const override {
        if (p.x < xvals_.front() - GEOM_EPS) {
            return INVALID_INDEX;
        }
        if (p.x > xvals_.back() + GEOM_EPS) {
            return INVALID_INDEX;
        }
        if (p.x < xvals_.front() + GEOM_EPS) {
            return 0;
        }

        return std::upper_bound(xvals_.begin(), xvals_.end(), p.x) - xvals_.begin() - 1;
    }

private:
    std::vector<double> xvals_;
};

//////////////////////////////////////////////////////////////////////
// 2D
//////////////////////////////////////////////////////////////////////
struct CellFinder::Impl::Finder2 : CellFinder::Impl {
    Finder2(const IGrid2D& g) : grid_(g), searcher_(build_boxes(g)) {}

    size_t find_cell(Point p) const override {
        for (size_t icell: searcher_.within(p)) {
            auto pp = grid_.tab_cell_point(icell) | std::views::transform([this](size_t i) { return grid_.point(i); });
            if (is_in_polygon(p, std::vector(pp.begin(), pp.end()), GEOM_EPS)) {
                return icell;
            }
        }
        return INVALID_INDEX;
    }

private:
    const IGrid2D& grid_;
    const BoxSearcher<2> searcher_;

    static std::vector<std::array<Point, 2>> build_boxes(const IGrid2D& g) {
        std::vector<std::array<Point, 2>> boxes;
        for (size_t icell = 0; icell < g.n_cells(); ++icell) {
            std::vector<size_t> ipoints = g.tab_cell_point(icell);
            std::array<Point, 2> box{g.point(ipoints[0]), g.point(ipoints[0])};
            for (size_t ipoint: ipoints | std::views::drop(1)) {
                Point p = g.point(ipoint);
                box[0].x = std::min(p.x, box[0].x);
                box[0].y = std::min(p.y, box[0].y);
                box[1].x = std::max(p.x, box[1].x);
                box[1].y = std::max(p.y, box[1].y);
            }
            boxes.push_back(box);
        }
        return boxes;
    }
};

CellFinder::CellFinder(const IGrid& grid) {
    switch (grid.dim()) {
    case 1:
        impl_ = std::make_unique<CellFinder::Impl::Finder1>(dynamic_cast<const IGrid1D&>(grid));
        break;
    case 2:
        impl_ = std::make_unique<CellFinder::Impl::Finder2>(dynamic_cast<const IGrid2D&>(grid));
        break;
    case 3:
        _THROW_NOT_IMP_;
    default:
        _THROW_UNREACHABLE_;
    };
}
CellFinder::~CellFinder() = default;

size_t CellFinder::operator()(Point p) const {
    return impl_->find_cell(p);
}
