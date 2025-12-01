#include "fem_assembler.hpp"
#include "cfd/mat/matrix_iter.hpp"

using namespace cfd;

FemAssembler::FemAssembler(size_t n_bases, const std::vector<FemElement>& elements,
                           const std::vector<std::vector<size_t>>& tab_elem_basis)
    : elements_(elements), tab_elem_basis_(tab_elem_basis) {

    // stencil
    std::vector<std::set<size_t>> tab_basis_basis(n_bases);
    for (size_t ielem = 0; ielem < elements.size(); ++ielem)
        for (size_t i = 0; i < tab_elem_basis[ielem].size(); ++i)
            for (size_t j = 0; j < i + 1; ++j) {
                size_t ibas1 = tab_elem_basis[ielem][i];
                size_t ibas2 = tab_elem_basis[ielem][j];
                tab_basis_basis[ibas1].insert(ibas2);
                tab_basis_basis[ibas2].insert(ibas1);
            }
    stencil_.set_stencil(tab_basis_basis);

    // element -> csr addresses table
    tab_elem_csr_address_.resize(n_elements());
    for (size_t ielem = 0; ielem < n_elements(); ++ielem) {
        size_t nbas = elements_[ielem].basis->size();
        for (size_t local_row = 0; local_row < nbas; ++local_row) {
            size_t global_row = tab_elem_basis[ielem][local_row];
            for (size_t local_col = 0; local_col < nbas; ++local_col) {
                size_t global_col = tab_elem_basis[ielem][local_col];
                size_t addr = stencil_.get_address(global_row, global_col);
                if (addr == INVALID_INDEX) {
                    _THROW_INTERNAL_ERROR_;
                }
                tab_elem_csr_address_[ielem].push_back(addr);
            }
        }
    }

    // basis types and reference points
    ref_points_.resize(n_bases);
    for (size_t ielem = 0; ielem < elements.size(); ++ielem) {
        auto geom = elements_[ielem].geometry;
        std::vector<Point> xi = elements_[ielem].basis->parametric_reference_points();
        for (size_t ibas = 0; ibas < xi.size(); ++ibas) {
            Point p = geom->to_physical(xi[ibas]);
            size_t ind = tab_elem_basis_[ielem][ibas];
            ref_points_[ind] = p;
        }
    }
}

size_t FemAssembler::n_elements() const {
    return elements_.size();
}

size_t FemAssembler::n_bases() const {
    return stencil_.n_rows();
}

const FemElement& FemAssembler::element(size_t ielem) const {
    return elements_[ielem];
}

const FemElement& FemAssembler::boundary_element(size_t iface) const {
    return boundary_elements_.at(iface);
}

Point FemAssembler::reference_point(size_t ibas) const {
    return ref_points_[ibas];
}

std::vector<size_t> FemAssembler::tab_elem_basis(size_t icell) const {
    return tab_elem_basis_[icell];
}

std::vector<double> FemAssembler::approximate(const IPointFunction& func) const {
    std::vector<double> ret;
    for (size_t ibas = 0; ibas < n_bases(); ++ibas) {
        Point p = reference_point(ibas);
        ret.push_back(func(p));
    }
    return ret;
}

std::vector<double> FemAssembler::local_vector(size_t ielem, const std::vector<double>& v) const {
    std::vector<double> ret;
    for (size_t bas: tab_elem_basis_[ielem]) {
        ret.push_back(v[bas]);
    }
    return ret;
}

std::vector<Vector> FemAssembler::local_vector(size_t ielem, const std::vector<double>& vx,
                                               const std::vector<double>& vy) const {

    std::vector<Vector> ret;
    for (size_t bas: tab_elem_basis_[ielem]) {
        ret.push_back({vx[bas], vy[bas], 0.0});
    }
    return ret;
}

std::vector<Vector> FemAssembler::local_vector(size_t ielem, const std::vector<double>& vx,
                                               const std::vector<double>& vy, const std::vector<double>& vz) const {

    std::vector<Vector> ret;
    for (size_t bas: tab_elem_basis_[ielem]) {
        ret.push_back({vx[bas], vy[bas], vz[bas]});
    }
    return ret;
}

std::vector<Vector> FemAssembler::local_vector(size_t ielem, const std::vector<Vector>& v) const {
    std::vector<Vector> ret;
    for (size_t bas: tab_elem_basis_[ielem]) {
        ret.push_back(v[bas]);
    }
    return ret;
}

void FemAssembler::add_to_global_matrix(double coef, size_t ielem, const std::vector<double>& local_matrix,
                                        std::vector<double>& global_csr_vals) const {
    for (size_t ival = 0; ival < local_matrix.size(); ++ival) {
        size_t a = tab_elem_csr_address_[ielem][ival];
        global_csr_vals[a] += coef * local_matrix[ival];
    }
}

void FemAssembler::add_to_global_matrix(size_t ielem, const std::vector<double>& local_matrix,
                                        std::vector<double>& global_csr_vals) const {
    add_to_global_matrix(1.0, ielem, local_matrix, global_csr_vals);
}

void FemAssembler::add_to_global_vector(double coef, size_t ielem, const std::vector<double>& local_vector,
                                        std::vector<double>& global_vector) const {
    for (size_t i = 0; i < local_vector.size(); ++i) {
        size_t gi = tab_elem_basis_[ielem][i];
        global_vector[gi] += coef * local_vector[i];
    }
}

void FemAssembler::add_to_global_vector(double coef, size_t ielem, const std::vector<double>& local_matrix,
                                        const std::vector<double>& u, std::vector<double>& global_vector) const {

    size_t n = tab_elem_basis_[ielem].size();
    for (size_t i = 0; i < n; ++i) {
        size_t gi = tab_elem_basis_[ielem][i];
        for (size_t j = 0; j < n; ++j) {
            size_t gj = tab_elem_basis_[ielem][j];
            global_vector[gi] += coef * local_matrix[n * i + j] * u[gj];
        }
    }
}

void FemAssembler::add_to_global_vector(size_t ielem, const std::vector<double>& local_vector,
                                        std::vector<double>& global_vector) const {
    add_to_global_vector(1, ielem, local_vector, global_vector);
}

const CsrStencil& FemAssembler::stencil() const {
    return stencil_;
}

void FemAssembler::boundary_add_to_global_vector(size_t iface, const std::vector<double>& local_vector,
                                                 std::vector<double>& global_vector) const {

    return add_to_global_vector(tab_face_elem_.at(iface), local_vector, global_vector);
}

void FemAssembler::add_boundary_element(size_t iface, size_t ielem, const FemElement& elem) {
    boundary_elements_[iface] = elem;
    tab_face_elem_[iface] = ielem;
}

CsrMatrix FemAssembler::zero_matrix() const {
    return CsrMatrix(stencil());
}

CsrMatrix FemAssembler::unit_matrix() const {
    CsrMatrix ret = zero_matrix();
    for (auto [d]: matrix_iter::v_diag(ret)) {
        d = 1.0;
    }
    return ret;
}
