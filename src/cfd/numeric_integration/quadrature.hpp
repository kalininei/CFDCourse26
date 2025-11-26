#ifndef __CFD_QUADRATURE_HPP__
#define __CFD_QUADRATURE_HPP__

#include "cfd/geom/point.hpp"
#include <array>
#include <functional>

namespace cfd {

class Quadrature {
public:
    Quadrature(const std::vector<Point>& points, const std::vector<double>& weights);

    size_t size() const;

    double integrate(const std::function<double(Point)>& func) const;
    double integrate(const std::vector<double>& values) const;

    std::vector<double> integrate(const std::function<std::vector<double>(Point)>& func) const;
    std::vector<double> integrate(const std::vector<std::vector<double>>& values) const;

    std::vector<Vector> integrate(const std::function<std::vector<Vector>(Point)>& func) const;
    std::vector<Vector> integrate(const std::vector<std::vector<Vector>>& values) const;

    const std::vector<Point>& points() const;
    const std::vector<double>& weights() const;

    Point point(size_t) const;
    double weight(size_t) const;

private:
    const std::vector<Point> points_;
    const std::vector<double> weights_;
};

}; // namespace cfd

#endif
