#ifndef __CFD_FEM_ASSEMBLER_HPP__
#define __CFD_FEM_ASSEMBLER_HPP__

#include "cfd/fem/fem_element.hpp"
#include "cfd/geom/i_point_function.hpp"
#include "cfd/grid/i_grid.hpp"
#include "cfd/mat/csrmat.hpp"

namespace cfd {

struct FemAssembler {
    FemAssembler(size_t n_bases, const std::vector<FemElement>& elements,
                 const std::vector<std::vector<size_t>>& tab_elem_basis);
    virtual ~FemAssembler() = default;

    void add_boundary_element(size_t iface, size_t ielem, const FemElement& elem);

    size_t n_elements() const;
    size_t n_bases() const;
    const CsrStencil& stencil() const;

    const FemElement& element(size_t icell) const;
    const FemElement& boundary_element(size_t iface) const;

    Point reference_point(size_t ibas) const;

    std::vector<size_t> tab_elem_basis(size_t icell) const;
    std::vector<double> approximate(const IPointFunction& func) const;
    std::vector<double> local_vector(size_t ielem, const std::vector<double>& v) const;
    std::vector<Vector> local_vector(size_t ielem, const std::vector<double>& vx, const std::vector<double>& vy) const;
    std::vector<Vector> local_vector(size_t ielem, const std::vector<double>& vx, const std::vector<double>& vy,
                                     const std::vector<double>& vz) const;

    void add_to_global_vector(size_t ielem, const std::vector<double>& local_vector,
                              std::vector<double>& global_vector) const;
    void add_to_global_vector(double coef, size_t ielem, const std::vector<double>& local_vector,
                              std::vector<double>& global_vector) const;
    void add_to_global_vector(double coef, size_t ielem, const std::vector<double>& local_matrix,
                              const std::vector<double>& u, std::vector<double>& global_vector) const;

    void add_to_global_matrix(size_t ielem, const std::vector<double>& local_matrix,
                              std::vector<double>& global_csr_vals) const;
    void add_to_global_matrix(double coef, size_t ielem, const std::vector<double>& local_matrix,
                              std::vector<double>& global_csr_vals) const;

    void boundary_add_to_global_vector(size_t iface, const std::vector<double>& local_vector,
                                       std::vector<double>& global_vector) const;

protected:
    std::vector<FemElement> _elements;

    std::vector<std::vector<size_t>> _tab_elem_basis;
    std::vector<std::vector<size_t>> _tab_elem_csr_address;
    std::vector<Point> _ref_points;
    std::map<size_t, size_t> tab_face_elem_;
    std::map<size_t, FemElement> boundary_elements_;

    CsrStencil _stencil;
};

} // namespace cfd
#endif
