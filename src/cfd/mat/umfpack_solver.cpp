#include "umfpack_solver.hpp"
#include <suitesparse/umfpack.h>
#include <stdexcept>
#include <cstring>
#include <iostream>

using namespace cfd;

struct UmfpackSolver::Impl {
    void* symbolic;
    void* numeric;
    
    std::vector<int> Ap;
    std::vector<int> Ai;
    std::vector<double> Ax;
    
    int n;
    bool matrix_initialized;
    bool factorized;
    
    Impl() : symbolic(nullptr), numeric(nullptr), n(0), 
             matrix_initialized(false), factorized(false) {}
    
    ~Impl() {
        if (numeric) {
            umfpack_di_free_numeric(&numeric);
        }
        if (symbolic) {
            umfpack_di_free_symbolic(&symbolic);
        }
    }
    
    void init_from_csr(const CsrMatrix& m) {
        if (numeric) {
            umfpack_di_free_numeric(&numeric);
            numeric = nullptr;
        }
        if (symbolic) {
            umfpack_di_free_symbolic(&symbolic);
            symbolic = nullptr;
        }
        
        n = static_cast<int>(m.n_rows());
        size_t nnz = m.n_nonzeros();
        
        if (n == 0 || nnz == 0) {
            throw std::runtime_error("Empty matrix provided");
        }

    // to CSC
    const auto& rowptr = m.addr();
    const auto& colind = m.cols();
    const auto& values = m.vals();

    // CSC формат: Ap (указатели колонок), Ai (индексы строк), Ax (значения)
    Ap.assign(n + 1, 0);  // размер n+1, заполняем нулями

    // Шаг 1: Подсчитываем количество элементов в каждой колонке
    for (size_t i = 0; i < nnz; ++i) {
        int col = static_cast<int>(colind[i]);
        Ap[col + 1]++;  // увеличиваем счетчик для следующей колонки
    }

    // Шаг 2: Вычисляем префиксные суммы (получаем начала колонок)
    for (int i = 0; i < n; ++i) {
        Ap[i + 1] += Ap[i];
    }

    // Шаг 3: Заполняем Ai и Ax
    Ai.resize(nnz);
    Ax.resize(nnz);

    // Вектор для отслеживания текущей позиции в каждой колонке
    std::vector<int> position = Ap;  // копия начал колонок

    // Проходим по всем ненулевым элементам исходной CSR матрицы
    for (int row = 0; row < n; ++row) {
        for (size_t j = rowptr[row]; j < rowptr[row + 1]; ++j) {
            int col = static_cast<int>(colind[j]);
            double val = values[j];

            // Вставляем в CSC формат
            int pos = position[col];     // текущая позиция в колонке col
            Ai[pos] = row;                // индекс строки
            Ax[pos] = val;                // значение
            position[col]++;               // двигаем указатель в этой колонке
        }
    }

    matrix_initialized = true;
    factorized = false;
}
        
    
    void factorize() {
        if (!matrix_initialized) {
            throw std::runtime_error("Matrix not set before factorize");
        }
        
        if (factorized) {
            return;
        }
        
        double control[UMFPACK_CONTROL];
        umfpack_di_defaults(control);
        
        control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_AUTO;
        control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
        
        int status = umfpack_di_symbolic(n, n, 
                                         Ap.data(), Ai.data(), Ax.data(),
                                         &symbolic, control, nullptr);
        
        if (status != UMFPACK_OK) {
            throw std::runtime_error("UMFPACK symbolic factorization failed: " + 
                                     std::to_string(status));
        }
        
        status = umfpack_di_numeric(Ap.data(), Ai.data(), Ax.data(),
                                    symbolic, &numeric, control, nullptr);
        
        if (status != UMFPACK_OK) {
            umfpack_di_free_symbolic(&symbolic);
            symbolic = nullptr;
            throw std::runtime_error("UMFPACK numeric factorization failed: " + 
                                     std::to_string(status));
        }
        
        factorized = true;
    }
    
    void solve_system(const std::vector<double>& rhs, std::vector<double>& x) {
        if (!matrix_initialized) {
            throw std::runtime_error("Matrix not set before solve");
        }
        
        if (!factorized) {
            throw std::runtime_error("Matrix is not factorized");
        }
        
        if (static_cast<int>(rhs.size()) != n) {
            throw std::runtime_error("RHS size does not match matrix dimension");
        }
        
        x.resize(n);
        
        int status = umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Ax.data(),
                                      x.data(), rhs.data(), numeric, nullptr, nullptr);
        
        if (status != UMFPACK_OK) {
            throw std::runtime_error("UMFPACK solve failed: " + std::to_string(status));
        }
    }
};

UmfpackSolver::UmfpackSolver() : pimpl_(std::make_unique<Impl>()) {}

UmfpackSolver::~UmfpackSolver() = default;

void UmfpackSolver::set_matrix(const CsrMatrix& m) {
    pimpl_->init_from_csr(m);
    pimpl_->factorize();
}

void UmfpackSolver::solve(const std::vector<double>& rhs, std::vector<double>& x) {
    pimpl_->solve_system(rhs, x);
}

void UmfpackSolver::solve_slae(const CsrMatrix& mat, const std::vector<double>& rhs, std::vector<double>& x){
    UmfpackSolver slv;
    slv.set_matrix(mat);
    slv.solve(rhs, x);
}
