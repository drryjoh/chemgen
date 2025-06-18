#include <iostream>
#include <array>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/IterativeSolvers>

int main() {
    constexpr int N = 5;
    using namespace Eigen;

    // std::array input
    std::array<std::array<double, N>, N> A_data = {{
        {{10, 2, 3, 4, 5}},
        {{1, 20, 3, 4, 5}},
        {{1, 2, 30, 4, 5}},
        {{1, 2, 3, 40, 5}},
        {{1, 2, 3, 4, 50}}
    }};
    std::array<double, N> b_data = {{130, 146, 177, 224, 290}};

    // SparseMatrix A
    SparseMatrix<double> A(N, N);

    // Fill using coeffRef or operator()(i,j)
    A.reserve(VectorXi::Constant(N, N));  // Optional: preallocate

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (A_data[i][j] != 0.0)
                A.coeffRef(i, j) = A_data[i][j];

    A.makeCompressed();  // Finalize structure

    // Fill b
    VectorXd b(N);
    for (int i = 0; i < N; ++i)
        b(i) = b_data[i];

    // GMRES solve
    GMRES<SparseMatrix<double>> solver;
    solver.setMaxIterations(50);
    solver.set_restart(5);
    solver.setTolerance(1e-12);
    solver.compute(A);
    VectorXd x = solver.solve(b);

    std::cout << "x = " << x.transpose() << "\n";
    std::cout << "Residual: " << (A * x - b).norm() << "\n";
    return 0;
}
