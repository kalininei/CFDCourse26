#include "quadrature.hpp"

using namespace cfd;

Quadrature::Quadrature(const std::vector<Point>& points, const std::vector<double>& weights)
    : points_(points), weights_(weights){};

const std::vector<Point>& Quadrature::points() const {
    return points_;
}

const std::vector<double>& Quadrature::weights() const {
    return weights_;
}

Point Quadrature::point(size_t i) const {
    return points_[i];
}
double Quadrature::weight(size_t i) const {
    return weights_[i];
}

size_t Quadrature::size() const {
    return points_.size();
}

double Quadrature::integrate(const std::vector<double>& values) const {
    double ret = 0;
    for (size_t i = 0; i < weights_.size(); ++i) {
        ret += weights_[i] * values[i];
    }
    return ret;
}

double Quadrature::integrate(const std::function<double(Point)>& func) const {
    std::vector<double> values;
    for (Point p : points_) {
        values.push_back(func(p));
    }
    return integrate(values);
}

std::vector<double> Quadrature::integrate(const std::vector<std::vector<double>>& values) const {
    const size_t n_out = values[0].size();
    std::vector<double> ret(n_out, 0);
    for (size_t i = 0; i < weights_.size(); ++i) {
        for (size_t j = 0; j < n_out; ++j) {
            ret[j] += weights_[i] * values[i][j];
        }
    }
    return ret;
}

std::vector<double> Quadrature::integrate(const std::function<std::vector<double>(Point)>& func) const {
    std::vector<std::vector<double>> values;
    for (Point p : points_) {
        values.push_back(func(p));
    }
    return integrate(values);
}
