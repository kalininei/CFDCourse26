#ifndef __CFD_FVM_DIVERGENCE_HPP__
#define __CFD_FVM_DIVERGENCE_HPP__

namespace cfd {

class FvmDivergence {
    FvmDivergence(const);

    // (ux, uy, uz) should be defined on collocation points
    std::vector<double> face_normals(const std::vector<double>& ux, const std::vector<double>& uy,
                                     const std::vector<double>& uz = {}) const;

private:
};

} // namespace cfd
#endif
