#include "fem_element.hpp"
#include <limits>

using namespace cfd;

struct FemElementValue::Cache {
    std::vector<double> phi;
    std::vector<Vector> grad_phi;
    std::vector<double> laplace;

    bool has_jac = false;
    JacobiMatrix jac;
};

FemElementValue::FemElementValue(const FemElement* elem, Point xi_point)
    : geom(elem->geometry.get()), basis(elem->basis.get()), xi(xi_point), pcache_(new Cache()) {}

FemElementValue::~FemElementValue(){};

double FemElementValue::phi(size_t i) const {
    if (pcache_->phi.size() == 0) {
        pcache_->phi = basis->value(xi);
    }
    return pcache_->phi[i];
}

Vector FemElementValue::grad_phi(size_t i) const {
    if (pcache_->grad_phi.size() == 0) {
        std::vector<Vector> grad_xi = basis->grad(xi);
        for (size_t j = 0; j < grad_xi.size(); ++j) {
            pcache_->grad_phi.push_back(gradient_to_physical(*jacobi(), grad_xi[j]));
        }
    }
    return pcache_->grad_phi[i];
}

double FemElementValue::laplace(size_t i) const {
    if (pcache_->laplace.size() == 0) {
        pcache_->laplace.resize(basis->size(), 0);
        double h = 1e-6; // <- TODO: Compute from jacobian

        Vector x = geom->to_physical(xi);
        Vector xi1_minus = geom->to_parametric(x - Vector{h, 0, 0});
        Vector xi1_plus = geom->to_parametric(x + Vector{h, 0, 0});
        Vector xi2_minus = geom->to_parametric(x - Vector{0, h, 0});
        Vector xi2_plus = geom->to_parametric(x + Vector{0, h, 0});
        Vector xi3_minus = geom->to_parametric(x - Vector{0, 0, h});
        Vector xi3_plus = geom->to_parametric(x + Vector{0, 0, h});

        std::vector<Vector> grad1_minus = basis->grad(xi1_minus);
        std::vector<Vector> grad1_plus = basis->grad(xi1_plus);
        std::vector<Vector> grad2_minus = basis->grad(xi2_minus);
        std::vector<Vector> grad2_plus = basis->grad(xi2_plus);
        std::vector<Vector> grad3_minus = basis->grad(xi3_minus);
        std::vector<Vector> grad3_plus = basis->grad(xi3_plus);

        for (size_t j = 0; j < basis->size(); ++j) {
            // use finite difference to compute divergence.
            // TODO: should be done using metric tensor and Hesse matrix
            double ddx0 = gradient_to_physical(*jacobi(), grad1_minus[j]).x;
            double ddx1 = gradient_to_physical(*jacobi(), grad1_plus[j]).x;
            double ddy0 = gradient_to_physical(*jacobi(), grad2_minus[j]).y;
            double ddy1 = gradient_to_physical(*jacobi(), grad2_plus[j]).y;
            double ddz0 = gradient_to_physical(*jacobi(), grad3_minus[j]).z;
            double ddz1 = gradient_to_physical(*jacobi(), grad3_plus[j]).z;
            pcache_->laplace[j] = (ddx1 - ddx0) / (2 * h) + (ddy1 - ddy0) / (2 * h) + (ddz1 - ddz0) / (2 * h);
        }
    }
    return pcache_->laplace[i];
}

const JacobiMatrix* FemElementValue::jacobi() const {
    if (!pcache_->has_jac) {
        pcache_->jac = geom->jacobi(xi);
        pcache_->has_jac = true;
    }
    return &(pcache_->jac);
}

double FemElementValue::modj() const {
    return jacobi()->modj;
}

double FemElementValue::interpolate(const std::vector<double>&) const {
    _THROW_NOT_IMP_;
}

Vector FemElementValue::interpolate(const std::vector<Vector>& f) const {
    Vector ret;
    for (size_t i = 0; i < basis->size(); ++i) {
        ret += phi(i) * f[i];
    }
    return ret;
}

double FemElementValue::divergence(const std::vector<Vector>& f) const {
    double ret = 0.0;
    for (size_t i = 0; i < basis->size(); ++i) {
        Vector dphi = grad_phi(i);
        ret += (f[i].x * dphi.x + f[i].y * dphi.y + f[i].z * dphi.z);
    }
    return ret;
}

namespace {

class BasedGeometry : public IElementGeometry {
public:
    BasedGeometry(std::shared_ptr<const IElementBasis> basis, const std::vector<Point>& x) : basis_(basis), x_(x) {
        if (basis_->size() != x_.size()) {
            throw std::runtime_error("Sizes do not match");
        }
    }

    JacobiMatrix jacobi(Point xi) const override {
        auto grad = basis_->grad(xi);

        Point d_dxi{};
        Point d_deta{};
        Point d_dzeta{};
        for (size_t i = 0; i < x_.size(); ++i) {
            d_dxi += grad[i].x * x_[i];
            d_deta += grad[i].y * x_[i];
            d_dzeta += grad[i].z * x_[i];
        };
        JacobiMatrix ret{
            d_dxi.x, d_deta.x, d_dzeta.x, //
            d_dxi.y, d_deta.y, d_dzeta.y, //
            d_dxi.z, d_deta.z, d_dzeta.z, //
        };
        fill_jacobi_modj_auto(ret);
        return ret;
    };

    Point to_physical(Point xi) const override {
        Point ret{};
        std::vector<double> v = basis_->value(xi);
        for (size_t i = 0; i < x_.size(); ++i) {
            ret += v[i] * x_[i];
        }

        return ret;
    }

private:
    std::shared_ptr<const IElementBasis> basis_;
    std::vector<Point> x_;
};
} // namespace

std::shared_ptr<IElementGeometry> cfd::build_geometry_from_basis(std::shared_ptr<IElementBasis> geometry_basis,
                                                                 const std::vector<Point>& x) {

    return std::make_shared<BasedGeometry>(geometry_basis, x);
}
