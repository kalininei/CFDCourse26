#include "cfd26_test.hpp"
#include "cfd/geom/searcher.hpp"

using namespace cfd;

TEST_CASE("Nearest point", "[nearest-point]"){
	std::vector<Point> points;
	points.push_back({0, 0, 0});
	points.push_back({1, 0, 0});
	points.push_back({0, 1, 0});
	points.push_back({1, 1, 0});
	points.push_back({0, 0, 1});

	PointSearcher<3> searcher(points);
	{
		std::vector<size_t> fnd = searcher.nearest({0, 0, 0}, 1);
		CHECK(fnd.size() == 1);
		CHECK(fnd[0] == 0);
	}
	{
		std::vector<size_t> fnd = searcher.nearest({0, 0.1, 0}, 2);
		std::sort(fnd.begin(), fnd.end());
		CHECK(fnd.size() == 2);
		CHECK(fnd[0] == 0);
		CHECK(fnd[1] == 2);
	}

	{
		std::vector<size_t> fnd = searcher.nearest({0, 0, 0.5}, 2);
		std::sort(fnd.begin(), fnd.end());
		CHECK(fnd.size() == 2);
		CHECK(fnd[0] == 0);
		CHECK(fnd[1] == 4);
	}
}
