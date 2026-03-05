#ifndef __CFD_MAT_MATRIX_ITER_HPP__
#define __CFD_MAT_MATRIX_ITER_HPP__

#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/i_mat.hpp"
#include <ranges>

namespace cfd {

namespace matrix_iter {

template<typename T>
concept CsrMatrixConcept = std::derived_from<std::remove_cvref_t<T>, CsrMatrix>;

// ========================= [ROW, COLUMN, VALUE]
// -> [irow, jcol, val&...]
template<CsrMatrixConcept M1, CsrMatrixConcept... M>
auto ijv(size_t start_index, size_t end_index, M1&& m1, M&&... m) {
    return std::views::iota(start_index, end_index) | std::views::transform([&](size_t irow) {
               return std::views::iota(m1.addr()[irow], m1.addr()[irow + 1]) |
                      std::views::transform([irow, &m1, &m...](size_t a) {
                          return std::make_tuple(std::cref(irow), std::cref(m1.cols()[a]), std::ref(m1.vals()[a]),
                                                 std::ref(m.vals()[a])...);
                      });
           }) |
           std::views::join;
}

template<CsrMatrixConcept M1, CsrMatrixConcept... M>
auto ijv(M1&& m1, M&&... m) {
    return ijv(0, m1.n_rows(), std::forward<M1>(m1), std::forward<M>(m)...);
}

// ========================== [COLUMN, VALUE]
// -> [jcol, val&...]
template<CsrMatrixConcept M1, CsrMatrixConcept... M>
auto jv(size_t irow, M1&& m1, M&&... m) {
    return std::views::iota(m1.addr()[irow], m1.addr()[irow + 1]) | std::views::transform([irow, &m1, &m...](size_t a) {
               return std::make_tuple(std::cref(m1.cols()[a]), std::ref(m1.vals()[a]), std::ref(m.vals()[a])...);
           });
}

// ========================== [VALUE]
// -> [val&...]
template<CsrMatrixConcept M1, CsrMatrixConcept... M>
auto v(M1&& m1, M&&... m) {
    return std::views::iota(0u, m1.n_nonzeros()) | std::views::transform([&](size_t a) {
               return std::make_tuple(std::ref(m1.vals()[a]), std::ref(m.vals()[a])...);
           });
};

// -> [val&...] for i == j
template<CsrMatrixConcept M1, CsrMatrixConcept... M>
auto v_diag(M1&& m1, M&&... m) {
    return std::views::iota(0u, m1.n_rows()) | std::views::transform([&](size_t irow) {
               const size_t a0 = m1.addr()[irow];
               const size_t a1 = m1.addr()[irow + 1];
               const auto b = m1.cols().cbegin();
               const size_t a = std::find(b + a0, b + a1, irow) - b;
               if (a == a1) {
                   throw std::runtime_error("no diagonal value in the matrix");
               };
               return std::make_tuple(std::ref(m1.vals()[a]), std::ref(m.vals()[a])...);
           });
};

// -> [irow, val&...] for i == j
template<CsrMatrixConcept M1, CsrMatrixConcept... M>
auto iv_diag(M1&& m1, M&&... m) {
    return std::views::iota(0u, m1.n_rows()) | std::views::transform([&](size_t irow) {
               const size_t a0 = m1.addr()[irow];
               const size_t a1 = m1.addr()[irow + 1];
               const auto b = m1.cols().cbegin();
               const size_t a = std::find(b + a0, b + a1, irow) - b;
               if (a == a1) {
                   throw std::runtime_error("no diagonal value in the matrix");
               };
               return std::make_tuple(irow, std::ref(m1.vals()[a]), std::ref(m.vals()[a])...);
           });
};

}; // namespace matrix_iter
}; // namespace cfd

#endif
