#include "fem_numeric_integrals.hpp"

using namespace cfd;

namespace {

std::vector<JacobiMatrix> build_quad_jacobi(const Quadrature& quad, const IElementGeometry& geom) {
    std::vector<JacobiMatrix> ret;

    for (Point xi : quad.points()) {
        ret.push_back(geom.jacobi(xi));
    }

    return ret;
}

std::vector<double> matrix_upper_to_sym(size_t nrows, const std::vector<double>& upper) {
    std::vector<double> ret(nrows * nrows);
    size_t k = 0;
    for (size_t irow = 0; irow < nrows; ++irow) {
        ret[irow * nrows + irow] = upper[k];
        k++;

        for (size_t icol = irow + 1; icol < nrows; ++icol) {
            ret[irow * nrows + icol] = upper[k];
            ret[icol * nrows + irow] = upper[k];
            k++;
        }
    }

    return ret;
}

} // namespace

NumericElementIntegrals::NumericElementIntegrals(const Quadrature* quad, std::shared_ptr<const IElementGeometry> geom,
                                                 std::shared_ptr<const IElementBasis> basis)
    : _quad(quad), _geom(geom), _basis(basis), _quad_jacobi(build_quad_jacobi(*quad, *geom)) {}

std::vector<double> NumericElementIntegrals::mass_matrix() const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<double> phi = _basis->value(quad_xi);

        std::vector<double> quad_mass;
        for (size_t irow = 0; irow < _basis->size(); ++irow)
            for (size_t icol = irow; icol < _basis->size(); ++icol) {
                quad_mass.push_back(phi[irow] * phi[icol] * _quad_jacobi[iquad].modj);
            }

        values.push_back(quad_mass);
    }

    std::vector<double> upper_mass = _quad->integrate(values);
    return matrix_upper_to_sym(_basis->size(), upper_mass);
}

std::vector<double> NumericElementIntegrals::load_vector() const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<double> phi = _basis->value(quad_xi);

        std::vector<double> quad_lvec;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            quad_lvec.push_back(phi[irow] * _quad_jacobi[iquad].modj);
        }

        values.push_back(quad_lvec);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::stiff_matrix() const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }

        std::vector<double> quad_stiff;
        for (size_t irow = 0; irow < _basis->size(); ++irow)
            for (size_t icol = irow; icol < _basis->size(); ++icol) {
                double d = dot_product(phi_grad_x[irow], phi_grad_x[icol]);
                quad_stiff.push_back(d * _quad_jacobi[iquad].modj);
            }

        values.push_back(quad_stiff);
    }

    std::vector<double> upper_stiff = _quad->integrate(values);
    return matrix_upper_to_sym(_basis->size(), upper_stiff);
}

std::vector<double> NumericElementIntegrals::transport_matrix(const std::vector<double>& vx,
                                                              const std::vector<double>& vy,
                                                              const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty())
                vel.x() += vx[i] * phi[i];
            if (!vy.empty())
                vel.y() += vy[i] * phi[i];
            if (!vz.empty())
                vel.z() += vz[i] * phi[i];
        }

        std::vector<double> quad_stiff;
        for (size_t irow = 0; irow < _basis->size(); ++irow)
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double d = dot_product(vel, phi_grad_x[icol]) * phi[irow];
                quad_stiff.push_back(d * _quad_jacobi[iquad].modj);
            }

        values.push_back(quad_stiff);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::mass_matrix_stab_supg(const std::vector<double>& vx,
                                                                   const std::vector<double>& vy,
                                                                   const std::vector<double>& vz) const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty())
                vel.x() += vx[i] * phi[i];
            if (!vy.empty())
                vel.y() += vy[i] * phi[i];
            if (!vz.empty())
                vel.z() += vz[i] * phi[i];
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow)
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double d = dot_product(vel, phi_grad_x[irow]) * phi[icol];
                q.push_back(d * _quad_jacobi[iquad].modj);
            }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::stiff_matrix_stab_supg(const std::vector<double>& vx,
                                                                    const std::vector<double>& vy,
                                                                    const std::vector<double>& vz) const {

    // TODO: now zero matrix for linear elements is returned

    size_t n = _basis->size();
    return std::vector<double>(n * n, 0.0);
}

std::vector<double> NumericElementIntegrals::transport_matrix_stab_supg(const std::vector<double>& vx,
                                                                        const std::vector<double>& vy,
                                                                        const std::vector<double>& vz) const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty())
                vel.x() += vx[i] * phi[i];
            if (!vy.empty())
                vel.y() += vy[i] * phi[i];
            if (!vz.empty())
                vel.z() += vz[i] * phi[i];
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double d1 = dot_product(vel, phi_grad_x[irow]);
            for (size_t icol = irow; icol < _basis->size(); ++icol) {
                double d2 = dot_product(vel, phi_grad_x[icol]);
                q.push_back(d1 * d2 * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    std::vector<double> upper = _quad->integrate(values);
    return matrix_upper_to_sym(_basis->size(), upper);
}

std::vector<double> NumericElementIntegrals::divergence_vector(const std::vector<double>& vx,
                                                               const std::vector<double>& vy,
                                                               const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        double div_u = 0;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty())
                div_u += vx[i] * phi_grad_x[i].x();
            if (!vy.empty())
                div_u += vy[i] * phi_grad_x[i].y();
            if (!vz.empty())
                div_u += vz[i] * phi_grad_x[i].z();
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            q.push_back(div_u * phi[irow] * _quad_jacobi[iquad].modj);
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::divergence_vector_byparts(const std::vector<double>& vx,
                                                                       const std::vector<double>& vy,
                                                                       const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel{0, 0, 0};
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty())
                vel.x() += vx[i] * phi[i];
            if (!vy.empty())
                vel.y() += vy[i] * phi[i];
            if (!vz.empty())
                vel.z() += vz[i] * phi[i];
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double v = vel.x() * phi_grad_x[irow].x() + vel.y() * phi_grad_x[irow].y() + vel.z() * phi_grad_x[irow].z();
            q.push_back(v * _quad_jacobi[iquad].modj);
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dx_matrix() const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].x() * phi[irow];
                q.push_back(v * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dy_matrix() const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].y() * phi[irow];
                q.push_back(v * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dz_matrix() const {
    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].z() * phi[irow];
                q.push_back(v * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dx_matrix_stab_supg2(const std::vector<double>& vx,
                                                                  const std::vector<double>& vy,
                                                                  const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel{0, 0, 0};
        ;
        double div_u = 0;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty()) {
                vel.x() += vx[i] * phi[i];
                div_u += vx[i] * phi_grad_x[i].x();
            }
            if (!vy.empty()) {
                vel.y() += vy[i] * phi[i];
                div_u += vy[i] * phi_grad_x[i].y();
            }
            if (!vz.empty()) {
                vel.z() += vz[i] * phi[i];
                div_u += vz[i] * phi_grad_x[i].z();
            }
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double p1 = dot_product(vel, phi_grad_x[irow]);
            double p2 = phi[irow] * div_u;
            double p = p1 + p2;
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].x();
                q.push_back(v * p * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dy_matrix_stab_supg2(const std::vector<double>& vx,
                                                                  const std::vector<double>& vy,
                                                                  const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel{0, 0, 0};
        ;
        double div_u = 0;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty()) {
                vel.x() += vx[i] * phi[i];
                div_u += vx[i] * phi_grad_x[i].x();
            }
            if (!vy.empty()) {
                vel.y() += vy[i] * phi[i];
                div_u += vy[i] * phi_grad_x[i].y();
            }
            if (!vz.empty()) {
                vel.z() += vz[i] * phi[i];
                div_u += vz[i] * phi_grad_x[i].z();
            }
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double p1 = dot_product(vel, phi_grad_x[irow]);
            double p2 = phi[irow] * div_u;
            double p = p1 + p2;
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].y();
                q.push_back(v * p * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dz_matrix_stab_supg2(const std::vector<double>& vx,
                                                                  const std::vector<double>& vy,
                                                                  const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel{0, 0, 0};
        ;
        double div_u = 0;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty()) {
                vel.x() += vx[i] * phi[i];
                div_u += vx[i] * phi_grad_x[i].x();
            }
            if (!vy.empty()) {
                vel.y() += vy[i] * phi[i];
                div_u += vy[i] * phi_grad_x[i].y();
            }
            if (!vz.empty()) {
                vel.z() += vz[i] * phi[i];
                div_u += vz[i] * phi_grad_x[i].z();
            }
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double p1 = dot_product(vel, phi_grad_x[irow]);
            double p2 = phi[irow] * div_u;
            double p = p1 + p2;
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].z();
                q.push_back(v * p * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dx_matrix_stab_supg(const std::vector<double>& vx,
                                                                 const std::vector<double>& vy,
                                                                 const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel{0, 0, 0};
        ;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty())
                vel.x() += vx[i] * phi[i];
            if (!vy.empty())
                vel.y() += vy[i] * phi[i];
            if (!vz.empty())
                vel.z() += vz[i] * phi[i];
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double p = dot_product(vel, phi_grad_x[irow]);
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].x();
                q.push_back(v * p * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dy_matrix_stab_supg(const std::vector<double>& vx,
                                                                 const std::vector<double>& vy,
                                                                 const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel{0, 0, 0};
        ;
        double div_u = 0;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty()) {
                vel.x() += vx[i] * phi[i];
                div_u += vx[i] * phi_grad_x[i].x();
            }
            if (!vy.empty()) {
                vel.y() += vy[i] * phi[i];
                div_u += vy[i] * phi_grad_x[i].y();
            }
            if (!vz.empty()) {
                vel.z() += vz[i] * phi[i];
                div_u += vz[i] * phi_grad_x[i].z();
            }
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double p = dot_product(vel, phi_grad_x[irow]);
            // double p2 = phi[irow] * div_u;
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].y();
                q.push_back(v * p * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::dz_matrix_stab_supg(const std::vector<double>& vx,
                                                                 const std::vector<double>& vy,
                                                                 const std::vector<double>& vz) const {

    std::vector<std::vector<double>> values;

    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<Vector> phi_grad_x;
        for (Vector phi_grad_xi : _basis->grad(quad_xi)) {
            phi_grad_x.push_back(gradient_to_physical(_quad_jacobi[iquad], phi_grad_xi));
        }
        std::vector<double> phi = _basis->value(quad_xi);

        Vector vel{0, 0, 0};
        double div_u = 0;
        for (size_t i = 0; i < _basis->size(); ++i) {
            if (!vx.empty()) {
                vel.x() += vx[i] * phi[i];
                div_u += vx[i] * phi_grad_x[i].x();
            }
            if (!vy.empty()) {
                vel.y() += vy[i] * phi[i];
                div_u += vy[i] * phi_grad_x[i].y();
            }
            if (!vz.empty()) {
                vel.z() += vz[i] * phi[i];
                div_u += vz[i] * phi_grad_x[i].z();
            }
        }

        std::vector<double> q;
        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            double p = dot_product(vel, phi_grad_x[irow]);
            // double p2 = phi[irow] * div_u;
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                double v = phi_grad_x[icol].z();
                q.push_back(v * p * _quad_jacobi[iquad].modj);
            }
        }

        values.push_back(q);
    }

    return _quad->integrate(values);
}

std::vector<double> NumericElementIntegrals::custom_matrix(OperandFunc f) const {
    std::vector<std::vector<double>> values;
    for (size_t iquad = 0; iquad < _quad->size(); ++iquad) {
        Point quad_xi = _quad->points()[iquad];
        std::vector<double> q;
        OperandArg arg(_geom.get(), _basis.get(), quad_xi);

        for (size_t irow = 0; irow < _basis->size(); ++irow) {
            for (size_t icol = 0; icol < _basis->size(); ++icol) {
                q.push_back(f(irow, icol, &arg));
            }
        }
        values.push_back(q);
    }
    return _quad->integrate(values);
}
