#ifndef __CFD_QUADRATURE_HPP__
#define __CFD_QUADRATURE_HPP__

#include <array>
#include <functional>
#include "cfd/geom/point.hpp"

namespace cfd{

class Quadrature{
public:
	Quadrature(const std::vector<Point>& points, const std::vector<double>& weights);

	size_t size() const;

	double integrate(const std::function<double(Point)>& func) const;
	double integrate(const std::vector<double>& values) const;

	std::vector<double> integrate(const std::function<std::vector<double>(Point)>& func) const;
	std::vector<double> integrate(const std::vector<std::vector<double>>& values) const;

	const std::vector<Point>& points() const;
	const std::vector<double>& weights() const;
private:
	const std::vector<Point> _points;
	const std::vector<double> _weights;
};

};

#endif
