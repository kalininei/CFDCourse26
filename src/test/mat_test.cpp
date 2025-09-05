#include "cfd26_test.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/sparse_matrix_solver.hpp"
#include "cfd/debug/printer.hpp"

using namespace cfd;

TEST_CASE("LodMatrix", "[lodmat]"){
	LodMatrix m(4);
	CHECK(m.n_rows() == 4);
	CHECK(m.n_nonzeros() == 0);

	// m[0, 1] = 0.12;
	m.add_value(0, 1, 0.12);
	CHECK(m.n_nonzeros() == 1);
	CHECK(m.value(0, 1) == 0.12);

	// m[0, 1] += 0.12;
	m.add_value(0, 1, 0.1);
	CHECK(m.n_nonzeros() == 1);
	CHECK(m.value(0, 1) == 0.22);

	//m[1, 1] = -0.5;
	m.add_value(1, 1, -0.5);
	CHECK(m.n_nonzeros() == 2);
	CHECK(m.value(1, 1) == -0.5);

	//m[0, :] = 0
	m.remove_row(0);
	CHECK(m.n_nonzeros() == 1);
	CHECK(m.value(0, 1) == 0);
	CHECK(m.value(1, 1) == -0.5);
};

TEST_CASE("CsrMatrix", "[csrmat]"){
	CsrMatrix m;

	// [1, *, *]
	// [*, 3, 1]
	// [1, *, 3]
	std::vector<size_t> addr{0, 1, 3, 5};
	std::vector<size_t> cols{0, 1, 2, 0, 2};
	std::vector<double> vals{1, 3, 1, 1, 3};

	m.set_stencil(std::move(addr), std::move(cols));
	m.set_values(std::move(vals));

	CHECK(m.n_rows() == 3);
	CHECK(m.n_nonzeros() == 5);
	CHECK(m.value(0, 0) == 1.0);
	CHECK(m.is_in_stencil(0, 0) == true);
	CHECK(m.is_in_stencil(0, 1) == false);
	CHECK(m.value(0, 1) == 0.0);
	CHECK(m.value(2, 2) == 3);

	// SLAE solution
	AmgcMatrixSolver solver;
	solver.set_matrix(m);
	std::vector<double> rhs{1, 1, 1};
	std::vector<double> x;
	solver.solve(rhs, x);

	CHECK(x[0] == Approx(1.0));
	CHECK(x[1] == Approx(0.333333));
	CHECK(x[2] == Approx(0));
}

TEST_CASE("Block matrix", "[block_matrix]"){

	// [1, *, *]
	// [*, 3, 1]
	// [1, *, 3]
	CsrMatrix M1(
		{0, 1, 3, 5},
		{0, 1, 2, 0, 2},
		{1, 3, 1, 1, 3}
	);

	{
		CsrMatrix r = assemble_block_matrix(3, 3, {{&M1}});
		CHECK(r.n_rows() == 3);
		CHECK(r.n_nonzeros() == 5);
		CHECK(r.value(0, 0) == 1);
		CHECK(r.value(2, 2) == 3);
		CHECK(r.value(1, 1) == 3);
		CHECK(r.value(0, 1) == 0);
	}
	{
		CsrMatrix r = assemble_block_matrix(3, 3, {{nullptr, &M1}});
		CHECK(r.n_rows() == 3);
		CHECK(r.n_nonzeros() == 5);
		CHECK(r.value(2, 2) == 0);
		CHECK(r.value(1, 1) == 0);
		CHECK(r.value(0, 1) == 0);
		CHECK(r.value(2, 5) == 3);
		CHECK(r.value(1, 4) == 3);
		CHECK(r.value(0, 4) == 0);
		CHECK(r.value(0, 3) == 1);
	}
	{
		CsrMatrix r = assemble_block_matrix(3, 3, {{nullptr}, {&M1}});
		CHECK(r.n_rows() == 6);
		CHECK(r.n_nonzeros() == 5);
		CHECK(r.value(2, 2) == 0);
		CHECK(r.value(1, 1) == 0);
		CHECK(r.value(0, 1) == 0);
		CHECK(r.value(5, 2) == 3);
		CHECK(r.value(4, 1) == 3);
		CHECK(r.value(3, 1) == 0);
		CHECK(r.value(3, 0) == 1);
	}
	{
		CsrMatrix r = assemble_block_matrix(3, 3, {{&M1, nullptr}, {&M1, &M1}});
		cfd::dbg::print(r);
		CHECK(r.n_rows() == 6);
		CHECK(r.n_nonzeros() == 15);
		CHECK(r.value(0, 0) == 1);
		CHECK(r.value(2, 2) == 3);
		CHECK(r.value(1, 1) == 3);
		CHECK(r.value(0, 1) == 0);
		CHECK(r.value(0, 3) == 0);
		CHECK(r.value(2, 5) == 0);
		CHECK(r.value(1, 4) == 0);
		CHECK(r.value(0, 4) == 0);
		CHECK(r.value(3, 0) == 1);
		CHECK(r.value(5, 2) == 3);
		CHECK(r.value(4, 1) == 3);
		CHECK(r.value(3, 1) == 0);
		CHECK(r.value(3, 3) == 1);
		CHECK(r.value(5, 5) == 3);
		CHECK(r.value(4, 4) == 3);
		CHECK(r.value(3, 4) == 0);
	}
}
