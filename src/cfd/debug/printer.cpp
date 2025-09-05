#include "printer.hpp"
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>

using namespace cfd;

void dbg::print(const ISparseMatrix& mat){
	size_t n = mat.n_rows();
	std::cout << "-- SIZE = " << n << "x" << n <<std::endl;
	for (size_t irow=0; irow < n; ++irow){
		for (size_t icol=0; icol < n; ++icol){
			std::cout << std::setw(6);
			if (mat.is_in_stencil(irow, icol)){
				std::cout << mat.value(irow, icol);
			} else {
				std::cout << "*";
			}
			std::cout << " ";
		}
		std::cout << std::endl;
	}
}

void dbg::print(size_t irow, const ISparseMatrix& mat){
	dbg::print(irow, mat, 0, mat.n_rows());
}

void dbg::print(size_t irow, const ISparseMatrix& mat, size_t col0, size_t col1){
	std::cout << "-- ROW = " << irow << std::endl;
	for (size_t icol=col0; icol < col1; ++icol){
		if (mat.is_in_stencil(irow, icol)){
			std::cout << "   [" << std::setw(4) << icol << "] ";
			std::cout << mat.value(irow, icol);
			std::cout << std::endl;
		}
	}
}

void dbg::print(const std::vector<double>& vec){
	std::cout << "-- SIZE = " << vec.size() << std::endl;
	for (double v: vec){
		std::cout << v << std::endl;
	}
}

void dbg::print_feat(const std::vector<double>& vec){
	double abs_sum = 0;
	double abs_min = std::numeric_limits<double>::max();
	for (double v: vec) {
		abs_sum += std::abs(v);
		abs_min = std::min(abs_min, std::abs(v));
	}
	std::cout << "SIZE:    " << vec.size() << std::endl;
	std::cout << "MIN:     " << *std::min_element(vec.begin(), vec.end()) << std::endl;
	std::cout << "MAX:     " << *std::max_element(vec.begin(), vec.end()) << std::endl;
	std::cout << "ABS_MIN: " << abs_min << std::endl;
	std::cout << "SUM:     " << std::accumulate(vec.begin(), vec.end(), 0.0) << std::endl;
	std::cout << "ABS_SUM: " << abs_sum << std::endl;
}

void dbg::ping_printer_cpp(){ }
