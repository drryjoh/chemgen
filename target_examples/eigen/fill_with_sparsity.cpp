//test
// File: sparse_gmres_fixed_pattern.cpp
#include <iostream>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/IterativeSolvers>

int main() {
    constexpr int N = 5;
    using namespace Eigen;

    // 1. Define known sparsity pattern (diagonal + first superdiagonal)
    SparseMatrix<double> A(N, N);
    A.reserve(VectorXi::Constant(N, 2));  // Each row has up to 2 entries

    for (int i = 0; i < N; ++i) {
        A.insert(i, i) = 0.0;       // diagonal
        if (i < N - 1) A.insert(i, i + 1) = 0.0;  // upper diagonal
    }

    A.makeCompressed();  // Finalize sparsity structure

    // 2. Define right-hand side (same for all iterations)
    VectorXd b(N);
    b << 1, 2, 3, 4, 5;

    // 3. Define GMRES solver
    GMRES<SparseMatrix<double>> solver;
    solver.setMaxIterations(20);
    solver.set_restart(5);
    solver.setTolerance(1e-10);

    // 4. Loop: refill A with new values and solve
    for (int iter = 0; iter < 5; ++iter) {
        // Overwrite all entries directly (no setZero() needed)
        for (int i = 0; i < N; ++i) {
            A.coeffRef(i, i) = 10 + iter + i;  // e.g., 10, 11, 12, ...
            if (i < N - 1) A.coeffRef(i, i + 1) = -1.0 * (iter + 1);
        }

        solver.compute(A);  // Pattern is reused, no overhead
        VectorXd x = solver.solve(b);

        // Output results
        std::cout << "Iteration " << iter << ": solution x = " << x.transpose() << "\n";
        std::cout << "Residual norm ||Ax - b|| = " << (A * x - b).norm() << "\n\n";
    }

    return 0;
}
