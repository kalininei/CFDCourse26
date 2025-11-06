#include "fvm_dfdn.hpp"
#include "cfd/debug/tictoc.hpp"
#include "cfd/geom/searcher.hpp"
#include "cfd/geom/simplex.hpp"
#include "cfd/mat/densemat.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// DuDn
///////////////////////////////////////////////////////////////////////////////

namespace {

LodMatrix assemble_faces_dudn_2d(const IGrid& grid, const FvmExtendedCollocations& colloc) {
    LodMatrix mat(grid.n_faces());

    std::vector<std::vector<size_t>> tab_point_colloc(grid.n_points());
    {
        for (size_t icolloc: colloc.cell_collocations) {
            size_t icell = colloc.cell_index(icolloc);
            for (size_t ipoint: grid.tab_cell_point(icell)) {
                tab_point_colloc[ipoint].push_back(icolloc);
            }
        }
        for (size_t icolloc: colloc.face_collocations) {
            size_t iface = colloc.face_index(icolloc);
            for (size_t ipoint: grid.tab_face_point(iface)) {
                tab_point_colloc[ipoint].push_back(icolloc);
            }
        }
    }
    auto find_closest_collocation = [&](size_t grid_point, size_t excl0, size_t excl1) {
        size_t ret = INVALID_INDEX;
        double min_meas = 1e100;
        Point p0 = grid.point(grid_point);
        for (size_t icolloc: tab_point_colloc[grid_point]) {
            if (icolloc != excl0 && icolloc != excl1) {
                double meas = vector_meas(p0 - colloc.points[icolloc]);
                if (meas < min_meas) {
                    min_meas = meas;
                    ret = icolloc;
                }
            }
        }
        return ret;
    };

    auto add_ds_entry = [&](const Vector& normal, const Vector& c, size_t col0, size_t col1, size_t col2,
                            size_t iface) {
        Vector s(-normal.y, normal.x);
        double cos_cos = dot_product(c, s) / dot_product(c, normal);

        Point p0 = colloc.points[col0];
        Point p1 = colloc.points[col1];
        Point p2 = colloc.points[col2];
        double tri_area = triangle_area(p0, p1, p2);

        double coef = -0.5 * cos_cos / tri_area;

        double x0 = p0.x;
        double y0 = p0.y;
        double x1 = p1.x;
        double y1 = p1.y;
        double x2 = p2.x;
        double y2 = p2.y;
        double dx0 = (y1 - y2) / 2.0;
        double dy0 = (x2 - x1) / 2.0;
        double dx1 = (y2 - y0) / 2.0;
        double dy1 = (x0 - x2) / 2.0;
        double dx2 = (y0 - y1) / 2.0;
        double dy2 = (x1 - x0) / 2.0;

        mat.add_value(iface, col0, coef * (dx0 * s.x + dy0 * s.y));
        mat.add_value(iface, col1, coef * (dx1 * s.x + dy1 * s.y));
        mat.add_value(iface, col2, coef * (dx2 * s.x + dy2 * s.y));
    };

    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        Vector normal = grid.face_normal(iface);
        size_t negative_collocation = colloc.tab_face_colloc(iface)[0];
        size_t positive_collocation = colloc.tab_face_colloc(iface)[1];
        Point ci = colloc.points[negative_collocation];
        Point cj = colloc.points[positive_collocation];

        // +dudc / cos(c,n);
        double v1 = 1.0 / dot_product(normal, cj - ci);
        mat.set_value(iface, positive_collocation, v1);
        mat.set_value(iface, negative_collocation, -v1);

        // -duds*cos(c,s)/cos(c,n)
        Vector c = (cj - ci) / vector_abs(cj - ci);
        {
            // left point
            size_t igrid = grid.tab_face_point(iface)[0];
            size_t col1 = find_closest_collocation(igrid, positive_collocation, negative_collocation);
            add_ds_entry(normal, c, negative_collocation, col1, positive_collocation, iface);
        }
        {
            // right point
            size_t igrid = grid.tab_face_point(iface)[1];
            size_t col1 = find_closest_collocation(igrid, positive_collocation, negative_collocation);
            add_ds_entry(normal, c, positive_collocation, col1, negative_collocation, iface);
        }
    }

    return mat;
}

LodMatrix build_dfdn_matrix(const IGrid& grid, const FvmExtendedCollocations& colloc) {
    if (grid.dim() == 2) {
        return assemble_faces_dudn_2d(grid, colloc);
    } else {
        _THROW_NOT_IMP_;
    }
}

} // namespace

FvmFacesDn::FvmFacesDn(const IGrid& grid) : FvmFacesDn(grid, FvmExtendedCollocations(grid)) {}

FvmFacesDn::FvmFacesDn(const IGrid& grid, const FvmExtendedCollocations& colloc)
    : dfdn_(build_dfdn_matrix(grid, colloc)) {}

std::vector<double> FvmFacesDn::compute(const std::vector<double>& f) const {
    return compute(f.data());
}

std::vector<double> FvmFacesDn::compute(const double* f) const {
    return dfdn_.mult_vec_p(f);
}

double FvmFacesDn::compute(size_t iface, const std::vector<double>& f) const {
    return compute(iface, f.data());
}

double FvmFacesDn::compute(size_t iface, const double* f) const {
    return dfdn_.mult_vec_p(iface, f);
}

const std::map<size_t, double>& FvmFacesDn::linear_combination(size_t iface) const {
    return dfdn_.row(iface);
}
