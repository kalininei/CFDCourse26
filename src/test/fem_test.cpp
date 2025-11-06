#include "cfd/fem/elem2d/quadrangle_linear.hpp"
#include "cfd/fem/fem_element.hpp"
#include "cfd26_test.hpp"

using namespace cfd;

TEST_CASE("Geometry from basis", "[geom-from-basis]") {
    {
        Point p0{0, 0};
        Point p1{1, 0.5};
        Point p2{1, 1.5};
        Point p3{0, 0.5};

        auto basis = std::make_shared<QuadrangleLinearBasis>();
        auto geom1 = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
        auto geom2 = build_geometry_from_basis(basis, {p0, p1, p2, p3});

        {
            double j1 = geom1->jacobi({}).modj;
            double j2 = geom2->jacobi({}).modj;
            CHECK(j1 == Approx(j2).margin(1e-8));
        }
        {
            double j1 = geom1->jacobi({0.1, 0.3}).modj;
            double j2 = geom2->jacobi({0.1, 0.3}).modj;
            CHECK(j1 == Approx(j2).margin(1e-8));
        }
        {
            double j1 = geom1->jacobi({1, -1}).modj;
            double j2 = geom2->jacobi({1, -1}).modj;
            CHECK(j1 == Approx(j2).margin(1e-8));
        }
    }

    {
        Point p0{0, 0};
        Point p1{1, 0.5};
        Point p2{0.8, 1.5};
        Point p3{0, 0.5};
        Point p4{0.5, 0.25};

        auto basis = std::make_shared<Quadrangle_Sqr0_Linear123>();
        auto geom1 = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
        auto geom2 = build_geometry_from_basis(basis, {p0, p1, p2, p3, p4});

        {
            double j1 = geom1->jacobi({}).modj;
            double j2 = geom2->jacobi({}).modj;
            CHECK(j1 == Approx(j2).margin(1e-8));
        }
        {
            double j1 = geom1->jacobi({0.1, 0.3}).modj;
            double j2 = geom2->jacobi({0.1, 0.3}).modj;
            CHECK(j1 == Approx(j2).margin(1e-8));
        }
        {
            double j1 = geom1->jacobi({1, -1}).modj;
            double j2 = geom2->jacobi({1, -1}).modj;
            CHECK(j1 == Approx(j2).margin(1e-8));
        }
    }

    {
        Point p0{0, 0};
        Point p1{1, 0.5};
        Point p2{1, 1.5};
        Point p3{0, 0.5};
        Point p4{0.5, 0.5};

        auto basis = std::make_shared<Quadrangle_Sqr0_Linear123>();
        auto geom1 = std::make_shared<QuadrangleLinearGeometry>(p0, p1, p2, p3);
        auto geom2 = build_geometry_from_basis(basis, {p0, p1, p2, p3, p4});

        Point xi0(0, -1);
        Point a0 = geom1->to_physical(xi0);
        Point b0 = geom2->to_physical(xi0);
        CHECK(vector_abs(a0 - Point(0.5, 0.25)) < 1e-8);
        CHECK(vector_abs(b0 - p4) < 1e-8);
    }
}
