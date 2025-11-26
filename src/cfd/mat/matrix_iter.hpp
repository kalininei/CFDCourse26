#ifndef __CFD_MAT_MATRIX_ITER_HPP__
#define __CFD_MAT_MATRIX_ITER_HPP__

#include "cfd/mat/csrmat.hpp"
#include "cfd/mat/i_mat.hpp"
#include <ranges>

namespace cfd {

namespace matrix_iter {

// ========================= [ROW, COLUMN, VALUE]
// -> [irow, jcol, val&]
inline auto ijv(CsrMatrix& m) {
    return std::views::iota(0u, m.n_rows()) | std::views::transform([&](size_t irow) {
               return std::views::iota(m.addr()[irow], m.addr()[irow + 1]) |
                      std::views::transform(
                          [irow, &m](size_t a) { return std::make_tuple(irow, m.cols()[a], std::ref(m.vals()[a])); });
           }) |
           std::views::join;
}
// -> [irow, jcol, val]
inline auto ijv(const CsrMatrix& m) {
    return std::views::iota(0u, m.n_rows()) | std::views::transform([&](size_t irow) {
               return std::views::iota(m.addr()[irow], m.addr()[irow + 1]) |
                      std::views::transform(
                          [irow, &m](size_t a) { return std::make_tuple(irow, m.cols()[a], m.vals()[a]); });
           }) |
           std::views::join;
}

// ========================== [COLUMN, VALUE]
// -> [jcol, val]
inline auto jv(size_t irow, const CsrMatrix& m) {
    return std::views::iota(m.addr()[irow], m.addr()[irow + 1]) |
           std::views::transform([&m](size_t a) { return std::make_tuple(m.cols()[a], m.vals()[a]); });
}
// -> [jcol, val&]
inline auto jv(size_t irow, CsrMatrix& m) {
    return std::views::iota(m.addr()[irow], m.addr()[irow + 1]) |
           std::views::transform([&m](size_t a) { return std::make_tuple(m.cols()[a], std::ref(m.vals()[a])); });
}
// -> [jcol, val1, val2]  !NOTE: m1 and m2 should have same stencil
inline auto jvv(size_t irow, const CsrMatrix& m1, const CsrMatrix& m2) {
    return std::views::iota(m1.addr()[irow], m1.addr()[irow + 1]) | std::views::transform([&m1, &m2](size_t a) {
               return std::make_tuple(m1.cols()[a], m1.vals()[a], m2.vals()[a]);
           });
}

}; // namespace matrix_iter
}; // namespace cfd

#endif
