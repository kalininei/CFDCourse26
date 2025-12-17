#include "cfd/grid/grid1d.hpp"
#include "cfd/grid/vtk.hpp"
#include "cfd26_test.hpp"

using namespace cfd;

namespace {

class BuckleyLeverettFunc {
public:
    double operator()(double x) {
        return x * x / (x * x + (1.0 - x) * (1.0 - x));
    };

    double dx(double x) {
        double d = x * x + (1.0 - x) * (1.0 - x);
        return -2 * x * (x - 1) / d / d;
    }
};

class LinearFunc {
public:
    LinearFunc(double x0, double x1, double v0, double v1) : x0_(x0), x1_(x1), v0_(v0), v1_(v1) {
        if (std::abs(x0 - x1) < 1e-12) {
            throw std::runtime_error("invalid linear function");
        }
    }

    double operator()(double x) {
        return (x - x0_) / (x1_ - x0_) * (v1_ - v0_) + v0_;
    }

    double dx(double) {
        return (v1_ - v0_) / (x1_ - x0_);
    }

private:
    const double x0_;
    const double x1_;
    const double v0_;
    const double v1_;
};

class ConcaveBuckleyLeverett {
public:
    ConcaveBuckleyLeverett() : f1_(), f2_(0, 0, lim_, f1_(lim_)){};

    double operator()(double x) {
        return (x < lim_) ? f2_(x) : f1_(x);
    }

    double dx(double x) {
        return (x < lim_) ? f2_.dx(x) : f1_.dx(x);
    }

private:
    const double lim_ = 0.7429341358783229;
    BuckleyLeverettFunc f1_;
    LinearFunc f2_;
};

class ConcaveBuckleyLeverett_Dx {
public:
    double operator()(double x) {
        return cbl_.dx(x);
    }

private:
    ConcaveBuckleyLeverett cbl_;
};

double inverse_monotonic_function(const ConcaveBuckleyLeverett_Dx&, double) {
    // TODO
    return 0;
};

class BuckleyLeverettSolver {
public:
    double solve(double x, double t) const {
        double xi = x / t;
        double a_left = tilde_f_->dx(u_left);
        double a_right = tilde_f_->dx(u_left);
        if (xi < a_left)
            return u_left;
        else if (xi > a_right)
            return u_right;
        else
            return inverse_monotonic_function(*tilde_dfdu_, xi);
    };

private:
    const double u_left = 1.0;
    const double u_right = 0.0;
    std::shared_ptr<ConcaveBuckleyLeverett> tilde_f_;
    std::shared_ptr<ConcaveBuckleyLeverett_Dx> tilde_dfdu_;
};

} // namespace

TEST_CASE("Buckley-Leverett exact solver", "[buckley-leverett-exact]") {
    BuckleyLeverettSolver slv;

    Grid1D g(0, 1, 100);

    VtkUtils::TimeSeriesWriter writer("bucklev");

    for (size_t i = 0; i < 10; ++i) {
        double time = 0.1 * static_cast<double>(i);
        std::string out_filename = writer.add(time);
        g.save_vtk(out_filename);

        std::vector<double> s(g.n_points(), 0.0);
        for (size_t ipoint = 0; ipoint < g.n_points(); ++ipoint) {
            s[i] = slv.solve(g.point(ipoint).x, time);
        }
        VtkUtils::add_point_data(s, "exact", out_filename);
    }
}
