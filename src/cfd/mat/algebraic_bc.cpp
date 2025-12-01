#include "algebraic_bc.hpp"
#include <ranges>

using namespace cfd;

// ================================ Dirichlet
void cfd::algebraic_bc_dirichlet(const std::map<size_t, double>& vals, CsrMatrix& lhs) {
    for (size_t i: vals | std::views::keys) {
        lhs.set_unit_row(i);
    }
}

void cfd::algebraic_bc_dirichlet(const std::map<size_t, double>& vals, std::vector<double>& rhs) {
    for (auto [i, v]: vals) {
        rhs[i] = v;
    }
}

// ================================ Periodic
void cfd::algebraic_bc_periodic(const std::map<size_t, size_t>& connections, CsrMatrix& lhs) {
    if (connections.empty())
        return;

    const size_t n_rows = lhs.addr().size() - 1;
    const auto& old_addr = lhs.addr();
    const auto& old_cols = lhs.cols();
    const auto& old_vals = lhs.vals();

    // 1. Affected rows
    std::vector<bool> is_affected(n_rows, false);
    for (const auto& [i, j]: connections) {
        is_affected[i] = true;
        is_affected[j] = true;
    }

    // 2. New stencil
    std::vector<size_t> new_row_sizes(n_rows, 0);
    std::vector<std::map<size_t, double>> sums_for_i;

    // not affected rows
    for (size_t row = 0; row < n_rows; row++) {
        if (!is_affected[row]) {
            new_row_sizes[row] = old_addr[row + 1] - old_addr[row];
        }
    }

    // affected rows
    for (auto [i, j]: connections) {
        // column->value for i-th rows
        std::map<size_t, double> sum_map;

        // i-th row
        size_t i_start = old_addr[i];
        size_t i_end = old_addr[i + 1];
        for (size_t idx = i_start; idx < i_end; idx++) {
            size_t col = old_cols[idx];
            double val = old_vals[idx];
            sum_map[col] = val;
        }

        // j-th row
        size_t j_start = old_addr[j];
        size_t j_end = old_addr[j + 1];
        for (size_t idx = j_start; idx < j_end; idx++) {
            size_t col = old_cols[idx];
            double val = old_vals[idx];
            if (col == j) {
                sum_map[i] += val; // A_ii += A_jj
            } else {
                sum_map[col] += val; // A_ik += A_jk
            }
        }

        new_row_sizes[i] = sum_map.size();
        new_row_sizes[j] = 2;
        sums_for_i.push_back(std::move(sum_map));
    }

    // 3. new csr
    std::vector<size_t> new_addr(n_rows + 1, 0);
    for (size_t row = 0; row < n_rows; row++) {
        new_addr[row + 1] = new_addr[row] + new_row_sizes[row];
    }
    size_t total_nnz = new_addr.back();
    std::vector<size_t> new_cols(total_nnz);
    std::vector<double> new_vals(total_nnz);

    // 4. copy not affected rows
    for (size_t row = 0; row < n_rows; row++) {
        if (is_affected[row])
            continue;

        size_t old_start = old_addr[row];
        size_t old_size = old_addr[row + 1] - old_start;
        size_t new_start = new_addr[row];

        if (old_size > 0) {
            std::copy(old_cols.begin() + old_start, old_cols.begin() + old_start + old_size,
                      new_cols.begin() + new_start);
            std::copy(old_vals.begin() + old_start, old_vals.begin() + old_start + old_size,
                      new_vals.begin() + new_start);
        }
    }

    // 5. copy affected rows
    size_t pair_index = 0;
    for (auto [i, j]: connections) {
        const auto& sum_map = sums_for_i[pair_index];

        // i-th rows
        size_t new_idx = new_addr[i];
        for (const auto& [col, val]: sum_map) {
            new_cols[new_idx] = col;
            new_vals[new_idx] = val;
            new_idx++;
        }

        // j-th rows
        size_t j_new_start1 = new_addr[j];
        size_t j_new_start2 = new_addr[j] + 1;
        if (j > i) {
            std::swap(j_new_start1, j_new_start2);
        }
        new_cols[j_new_start1] = j;
        new_vals[j_new_start1] = 1.0;
        new_cols[j_new_start2] = i;
        new_vals[j_new_start2] = -1.0;

        pair_index++;
    }

    lhs = CsrMatrix(std::move(new_addr), std::move(new_cols), std::move(new_vals));
}

void cfd::algebraic_bc_periodic(const std::map<size_t, size_t>& vals, std::vector<double>& rhs) {
    for (auto [i, j]: vals) {
        rhs[i] += rhs[j];
        rhs[j] = 0.0;
    }
}
