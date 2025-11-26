#include "cfd/geom/polygon.hpp"
#include "cfd/geom/searcher.hpp"
#include "cfd26_test.hpp"

using namespace cfd;

TEST_CASE("Nearest point", "[nearest-point]") {
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

TEST_CASE("Point in polygon with tolerance", "[is_in_polygon]") {
    std::vector<Point> triangle_ccw = {{0, 0}, {1, 0}, {0, 1}};
    std::vector<Point> square_ccw = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

    SECTION("Points strictly inside") {
        REQUIRE(is_in_polygon({0.3, 0.3}, triangle_ccw));
        REQUIRE(is_in_polygon({0.1, 0.1}, triangle_ccw));
        REQUIRE(is_in_polygon({0.5, 0.5}, square_ccw));
        REQUIRE(is_in_polygon({0.1, 0.9}, square_ccw));
    }

    SECTION("Points strictly outside") {
        REQUIRE_FALSE(is_in_polygon({0.8, 0.8}, triangle_ccw));
        REQUIRE_FALSE(is_in_polygon({-0.1, 0.5}, triangle_ccw));
        REQUIRE_FALSE(is_in_polygon({1.5, 0.5}, square_ccw));
        REQUIRE_FALSE(is_in_polygon({0.5, -0.1}, square_ccw));
    }

    SECTION("Points on boundary") {
        REQUIRE(is_in_polygon({0.5, 0.0}, triangle_ccw));
        REQUIRE(is_in_polygon({0.0, 0.5}, triangle_ccw));
        REQUIRE(is_in_polygon({0.5, 0.5}, triangle_ccw));
        REQUIRE(is_in_polygon({0.5, 0.0}, square_ccw));
        REQUIRE(is_in_polygon({1.0, 0.5}, square_ccw));
    }

    SECTION("Points near boundary within tolerance") {
        double eps = 1e-12;
        REQUIRE(is_in_polygon({0.5 + eps / 2, 0.0}, triangle_ccw, eps));
        REQUIRE(is_in_polygon({0.0, 0.5 + eps / 2}, triangle_ccw, eps));
        REQUIRE(is_in_polygon({0.5, -eps / 2}, square_ccw, eps));
        REQUIRE(is_in_polygon({1.0 + eps / 2, 0.5}, square_ccw, eps));
    }

    SECTION("Points near boundary outside tolerance") {
        double eps = 1e-12;
        REQUIRE_FALSE(is_in_polygon({1.0 + eps * 2, 0.0}, triangle_ccw, eps));
        REQUIRE_FALSE(is_in_polygon({0.0, -eps * 2}, square_ccw, eps));
    }

    SECTION("Vertex points") {
        REQUIRE(is_in_polygon({0, 0}, triangle_ccw));
        REQUIRE(is_in_polygon({1, 0}, triangle_ccw));
        REQUIRE(is_in_polygon({0, 1}, triangle_ccw));
        REQUIRE(is_in_polygon({0, 0}, square_ccw));
        REQUIRE(is_in_polygon({1, 1}, square_ccw));
    }

    SECTION("Custom tolerance") {
        Point near_boundary = {0.5, -0.5e-12};

        REQUIRE(is_in_polygon(near_boundary, square_ccw, 1e-12));

        REQUIRE_FALSE(is_in_polygon(near_boundary, square_ccw, 1e-13));

        REQUIRE(is_in_polygon({0.5, -1e-10}, square_ccw, 1e-9));
    }

    SECTION("Complex polygon") {
        std::vector<Point> hexagon = {{0, 0}, {2, 0}, {3, 1}, {2, 2}, {0, 2}, {-1, 1}};

        REQUIRE(is_in_polygon({1, 1}, hexagon));
        REQUIRE(is_in_polygon({2, 1}, hexagon));
        REQUIRE_FALSE(is_in_polygon({3, 2}, hexagon));
        REQUIRE_FALSE(is_in_polygon({0, 3}, hexagon));
    }

    SECTION("Edge cases") {
        std::vector<Point> empty_polygon;
        REQUIRE_FALSE(is_in_polygon({0, 0}, empty_polygon));

        std::vector<Point> line = {{0, 0}, {1, 0}};
        REQUIRE_FALSE(is_in_polygon({0.5, 0.0}, line));
        REQUIRE_FALSE(is_in_polygon({0.5, 0.1}, line));

        REQUIRE_FALSE(is_in_polygon({100, 100}, triangle_ccw));
        REQUIRE_FALSE(is_in_polygon({-100, -100}, triangle_ccw));
    }

    SECTION("Zero tolerance") {
        REQUIRE(is_in_polygon({0.5, 1e-13}, triangle_ccw, 0.0));
        REQUIRE_FALSE(is_in_polygon({0.5, -1e-13}, triangle_ccw, 0.0));
    }

    SECTION("Large tolerance") {
        REQUIRE(is_in_polygon({0.5, -0.1}, square_ccw, 0.2));
        REQUIRE(is_in_polygon({1.1, 0.5}, square_ccw, 0.2));
    }
}

TEST_CASE("Point within bbox", "[bbox-searcher]") {
    std::vector<std::array<Point, 2>> boxes;
    boxes.push_back(std::array<Point, 2>{Point{0.0, 0.0}, Point{1.0, 1.0}});
    boxes.push_back(std::array<Point, 2>{Point{0.5, 0.5}, Point{1.5, 1.5}});
    boxes.push_back(std::array<Point, 2>{Point{-0.5, 0.5}, Point{0.5, 1.5}});
    boxes.push_back(std::array<Point, 2>{Point{1.0, 1.0}, Point{2.0, 2.0}});

    BoxSearcher<2> searcher(boxes);
    {
        std::vector<size_t> fnd = searcher.within({0, 0});
        CHECK(fnd == std::vector<size_t>{0});
    }
    {
        std::vector<size_t> fnd = searcher.within({0.1, 0.1});
        CHECK(fnd == std::vector<size_t>{0});
    }
    {
        std::vector<size_t> fnd = searcher.within({0.99, 0.99});
        CHECK(fnd == std::vector<size_t>{0, 1});
    }
    {
        std::vector<size_t> fnd = searcher.within({1.0, 1.0});
        CHECK(fnd == std::vector<size_t>{0, 1, 3});
    }
    {
        std::vector<size_t> fnd = searcher.within({0.5, 0.7});
        CHECK(fnd == std::vector<size_t>{0, 1, 2});
    }
}
