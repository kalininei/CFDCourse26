#include "vecmat.hpp"

std::vector<double> compute_residual_vec(
		const cfd::CsrMatrix& m,
		const std::vector<double>& rhs,
		const std::vector<double>& u){

	std::vector<double> mvec = m.mult_vec(u);
	for (size_t i=0; i<mvec.size(); ++i){
		mvec[i] = std::abs(mvec[i] - rhs[i]);
	}
	return mvec;
}

double compute_residual(const cfd::CsrMatrix& m, const std::vector<double>& rhs, const std::vector<double>& u){
	std::vector<double> res = compute_residual_vec(m, rhs, u);
	return *std::max_element(res.begin(), res.end());
}

std::vector<double> vector_sum(
		const std::vector<double>& v1,
		double coef,
		const std::vector<double>& v2){

	std::vector<double> ret(v1);
	for (size_t i=0; i<ret.size(); ++i){
		ret[i] += coef * v2[i];
	}
	return ret;
}

std::vector<cfd::Vector> vector_sum(
		const std::vector<cfd::Vector>& v1,
		double coef,
		const std::vector<cfd::Vector>& v2){

	std::vector<cfd::Vector> ret(v1);
	for (size_t i=0; i<ret.size(); ++i){
		ret[i] += coef * v2[i];
	}
	return ret;
}

double harmonic_mean(double a, double b){
	return 2*a*b/(a+b);
}
