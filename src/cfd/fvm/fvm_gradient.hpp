#ifndef CFD_FVM_GRADIENT_HPP
#define CFD_FVM_GRADIENT_HPP

#include "cfd/fvm/fvm_extended_collocations.hpp"
#include "cfd/mat/lodmat.hpp"

namespace cfd {

struct IFvmGradient {
    virtual ~IFvmGradient() = default;

    std::vector<Vector> compute(const std::vector<double>& u) const;
    std::vector<Vector> compute(const double* u) const;

    const CsrMatrix& mat_x() const {
        return data_[0];
    }
    const CsrMatrix& mat_y() const {
        return data_[1];
    }
    const CsrMatrix& mat_z() const {
        return data_[2];
    }

protected:
    std::array<CsrMatrix, 3> data_;
};

///////////////////////////////////////////////////////////////////////////////
// FvmCellGradient
///////////////////////////////////////////////////////////////////////////////

struct IFvmCellGradient : public IFvmGradient {
    virtual ~IFvmCellGradient() = default;
};

struct LeastSquaresFvmCellGradient : public IFvmCellGradient {
    LeastSquaresFvmCellGradient(const IGrid& grid, const FvmExtendedCollocations& colloc);
};

struct GaussLinearFvmCellGradient : public IFvmCellGradient {
    GaussLinearFvmCellGradient(const IGrid& grid, const FvmExtendedCollocations& colloc);
};

///////////////////////////////////////////////////////////////////////////////
// FvmFaceGradient
///////////////////////////////////////////////////////////////////////////////

struct IFvmFaceGradient : public IFvmGradient {
    virtual ~IFvmFaceGradient() = default;

protected:
    IFvmFaceGradient(const IGrid& grid, const FvmExtendedCollocations& colloc, const IFvmCellGradient& cg);
};

struct LeastSquaresFvmFaceGradient : public IFvmFaceGradient {
    LeastSquaresFvmFaceGradient(const IGrid& grid, const FvmExtendedCollocations& colloc)
        : IFvmFaceGradient(grid, colloc, LeastSquaresFvmCellGradient(grid, colloc)) {}
};

struct GaussLinearFvmFaceGradient : public IFvmFaceGradient {
    GaussLinearFvmFaceGradient(const IGrid& grid, const FvmExtendedCollocations& colloc)
        : IFvmFaceGradient(grid, colloc, GaussLinearFvmCellGradient(grid, colloc)) {}
};

} // namespace cfd
#endif
