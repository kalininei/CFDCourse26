#include "quadrature.hpp"

using namespace cfd;

Quadrature::Quadrature(const std::vector<Point>& points, const std::vector<double>& weights)
	:_points(points), _weights(weights) { };

const std::vector<Point>& Quadrature::points() const{
	return _points;
}

const std::vector<double>& Quadrature::weights() const{
	return _weights;
}

size_t Quadrature::size() const{
	return _points.size();
}

double Quadrature::integrate(const std::vector<double>& values) const {
	double ret = 0;
	for (size_t i=0; i<_weights.size(); ++i){
		ret += _weights[i] * values[i];
	}
	return ret;
}

double Quadrature::integrate(const std::function<double(Point)>& func) const{
	std::vector<double> values;
	for (Point p: _points){
		values.push_back(func(p));
	}
	return integrate(values);
}

std::vector<double> Quadrature::integrate(const std::vector<std::vector<double>>& values) const {
	const size_t n_out = values[0].size();
	std::vector<double> ret(n_out, 0);
	for (size_t i=0; i<_weights.size(); ++i){
		for (size_t j=0; j<n_out; ++j){
			ret[j] += _weights[i] * values[i][j];
		}
	}
	return ret;
}

std::vector<double> Quadrature::integrate(const std::function<std::vector<double>(Point)>& func) const{
	std::vector<std::vector<double>> values;
	for (Point p: _points){
		values.push_back(func(p));
	}
	return integrate(values);
}
