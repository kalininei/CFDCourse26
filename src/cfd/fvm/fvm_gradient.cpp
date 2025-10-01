#include "fvm_gradient.hpp"
#include "cfd/debug/tictoc.hpp"
#include "cfd/mat/densemat.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Cell gradient
///////////////////////////////////////////////////////////////////////////////

namespace {

DenseMatrix least_squares_inv(const DenseMatrix& a) {
    // transpose(A)
    DenseMatrix at = a.transpose();
    // inverse(transpose(A) * A)
    DenseMatrix inv_at_a = at.mult_mat(a).inverse();
    // inverse(transpose(A) * A) * transpose(A)
    return inv_at_a.mult_mat(at);
}

std::array<CsrMatrix, 3> assemble_fvm_cell_gradient_2d(const IGrid& grid, const FvmExtendedCollocations& colloc) {
    LodMatrix grad_x(grid.n_cells());
    LodMatrix grad_y(grid.n_cells());

    for (size_t icell = 0; icell < grid.n_cells(); ++icell) {
        const std::vector<size_t>& collocs = colloc.tab_colloc_colloc(icell);

        DenseMatrix amat(collocs.size(), 2);
        for (size_t i = 0; i < collocs.size(); ++i) {
            Vector c = colloc.points[collocs[i]] - colloc.points[icell];
            amat.set_value(i, 0, c.x);
            amat.set_value(i, 1, c.y);
        }
        DenseMatrix lsi = least_squares_inv(amat);
        double diag_x = 0;
        double diag_y = 0;
        for (size_t i = 0; i < collocs.size(); ++i) {
            double vx = lsi.value(0, i);
            double vy = lsi.value(1, i);
            grad_x.set_value(icell, collocs[i], vx);
            grad_y.set_value(icell, collocs[i], vy);
            diag_x -= vx;
            diag_y -= vy;
        }
        grad_x.set_value(icell, icell, diag_x);
        grad_y.set_value(icell, icell, diag_y);
    }

    return {grad_x.to_csr(), grad_y.to_csr()};
}

std::array<CsrMatrix, 3> assemble_fvm_cell_gradient(const IGrid& grid, const FvmExtendedCollocations& colloc) {
    if (grid.dim() == 2) {
        return assemble_fvm_cell_gradient_2d(grid, colloc);
    } else {
        _THROW_NOT_IMP_;
    }
}

} // namespace

LeastSquaresFvmCellGradient::LeastSquaresFvmCellGradient(const IGrid& grid, const FvmExtendedCollocations& colloc) {
    data_ = assemble_fvm_cell_gradient(grid, colloc);
}

GaussLinearFvmCellGradient::GaussLinearFvmCellGradient(const IGrid& grid, const FvmExtendedCollocations& colloc) {
    LodMatrix grad_x(grid.n_cells());
    LodMatrix grad_y(grid.n_cells());
    LodMatrix grad_z(grid.n_cells());

    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        auto icollocs = colloc.tab_face_colloc(iface);
        Vector n = grid.face_area(iface) * grid.face_normal(iface);
        // ---- internal face
        if (icollocs[0] < grid.n_cells() && icollocs[1] < grid.n_cells()) {
            const auto& icells = icollocs;
            double vol1 = grid.cell_volume(icells[0]);
            double vol2 = grid.cell_volume(icells[1]);
            // left cell
            grad_x.add_value(icells[0], icells[0], n.x * 0.5 / vol1);
            grad_y.add_value(icells[0], icells[0], n.y * 0.5 / vol1);
            grad_z.add_value(icells[0], icells[0], n.z * 0.5 / vol1);
            grad_x.add_value(icells[0], icells[1], n.x * 0.5 / vol1);
            grad_y.add_value(icells[0], icells[1], n.y * 0.5 / vol1);
            grad_z.add_value(icells[0], icells[1], n.z * 0.5 / vol1);
            // right cell
            grad_x.add_value(icells[1], icells[0], -n.x * 0.5 / vol2);
            grad_y.add_value(icells[1], icells[0], -n.y * 0.5 / vol2);
            grad_z.add_value(icells[1], icells[0], -n.z * 0.5 / vol2);
            grad_x.add_value(icells[1], icells[1], -n.x * 0.5 / vol2);
            grad_y.add_value(icells[1], icells[1], -n.y * 0.5 / vol2);
            grad_z.add_value(icells[1], icells[1], -n.z * 0.5 / vol2);
        }
        // boundary face (no right cell)
        else if (icollocs[1] >= grid.n_cells()) {
            double vol1 = grid.cell_volume(icollocs[0]);
            grad_x.add_value(icollocs[0], icollocs[1], n.x / vol1);
            grad_y.add_value(icollocs[0], icollocs[1], n.y / vol1);
            grad_z.add_value(icollocs[0], icollocs[1], n.z / vol1);
            // boundary face (no left cell)
        } else {
            double vol2 = grid.cell_volume(icollocs[1]);
            grad_x.add_value(icollocs[1], icollocs[0], -n.x / vol2);
            grad_y.add_value(icollocs[1], icollocs[0], -n.y / vol2);
            grad_z.add_value(icollocs[1], icollocs[0], -n.z / vol2);
        }
    }
    data_[0] = grad_x.to_csr();
    data_[1] = grad_y.to_csr();
    if (grid.dim() > 2) {
        data_[2] = grad_z.to_csr();
    }
}

std::vector<Vector> IFvmGradient::compute(const double* u) const {
    std::vector<double> x = data_[0].mult_vec_p(u);
    std::vector<double> y = data_[1].mult_vec_p(u);

    std::vector<Vector> ret(x.size());
    for (size_t i = 0; i < ret.size(); ++i) {
        ret[i].x = x[i];
        ret[i].y = y[i];
        ret[i].z = 0.0;
    }

    if (data_[2].n_rows() > 0) {
        std::vector<double> z = data_[2].mult_vec_p(u);
        for (size_t i = 0; i < ret.size(); ++i) {
            ret[i].z = z[i];
        }
    }

    return ret;
}

std::vector<Vector> IFvmGradient::compute(const std::vector<double>& u) const {
    return compute(u.data());
}

IFvmFaceGradient::IFvmFaceGradient(const IGrid& grid, const FvmExtendedCollocations& colloc,
                                   const IFvmCellGradient& cg) {
    LodMatrix grad_x(grid.n_faces());
    LodMatrix grad_y(grid.n_faces());
    LodMatrix grad_z(grid.n_faces());

    auto copy_row = [&](double c, size_t target, size_t irow) {
        size_t astart = cg.mat_x().addr()[irow];
        size_t aend = cg.mat_x().addr()[irow + 1];
        const size_t* cols = cg.mat_x().cols().data() + astart;
        const double* xdata = cg.mat_x().vals().data() + astart;
        const double* ydata = cg.mat_y().vals().data() + astart;
        const double* zdata = (grid.dim() > 2) ? cg.mat_z().vals().data() + astart : nullptr;
        for (size_t a = 0; a < aend - astart; ++a) {
            size_t col = cols[a];
            grad_x.add_value(target, col, c * xdata[a]);
            grad_y.add_value(target, col, c * ydata[a]);
            if (grid.dim() > 2) {
                grad_z.add_value(target, col, c * zdata[a]);
            }
        }
    };

    for (size_t iface = 0; iface < grid.n_faces(); ++iface) {
        auto [i, j] = colloc.tab_face_colloc(iface);
        if (colloc.is_boundary_colloc(i)) {
            copy_row(1.0, iface, j);
        } else if (colloc.is_boundary_colloc(j)) {
            copy_row(1.0, iface, i);
        } else {
            double wi = 0.5;
            double wj = 0.5;
            copy_row(wi, iface, i);
            copy_row(wj, iface, j);
        }
    }

    data_[0] = grad_x.to_csr();
    data_[1] = grad_y.to_csr();
    if (grid.dim() > 2) {
        data_[2] = grad_z.to_csr();
    }
}
