//Generated
// File: gmres_example.cpp
#include <iostream>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/IterativeSolvers>

int main() {
    using namespace Eigen;

    // Define a sparse matrix A
    typedef SparseMatrix<double> SpMat;
    typedef Triplet<double> T;
    std::vector<T> triplets;

    const int size = 5;
    for (int i = 0; i < size; ++i) {
        triplets.emplace_back(i, i, 4.0); // diagonal
        if (i > 0) triplets.emplace_back(i, i - 1, -1.0); // lower diagonal
        if (i < size - 1) triplets.emplace_back(i, i + 1, -1.0); // upper diagonal
    }

    SpMat A(size, size);
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Right-hand side vector b
    VectorXd b(size);
    b << 1, 2, 3, 4, 5;

    // Initial guess (zero)
    VectorXd x = VectorXd::Zero(size);

    // GMRES solver (Eigen's unsupported module)
    Eigen::GMRES<SpMat> solver;
    solver.setMaxIterations(100);
    solver.set_restart(10);  // GMRES restart value
    solver.setTolerance(1e-10);
    solver.compute(A);

    x = solver.solve(b);

    std::cout << "Solution x:\n" << x << "\n";
    std::cout << "Number of iterations: " << solver.iterations() << "\n";
    std::cout << "Estimated error: " << solver.error() << "\n";

    return 0;
}
