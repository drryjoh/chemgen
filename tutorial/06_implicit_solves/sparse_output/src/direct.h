#ifdef CHEMGEN_EIGEN


ChemicalState
direct_solve(const SparseMatrix<double>& A, const ChemicalState& b_, SparseLU<SparseMatrix<double>, COLAMDOrdering<int>>& lu, bool compute)
{
    Matrix<double, n_variables, 1> b(n_variables);

    for (int i = 0; i < n_variables; ++i)
    {
        b(i) = b_[i];
    }

    Matrix<double, n_variables, 1> x;

#ifdef CHEMGEN_EIGEN_SPARSE

    if (compute) lu.compute(A);
    if (lu.info() != Success) std::cerr << "lu factorization failed!" << std::endl;
    x = lu.solve(b);

#else

    if (compute) lu = SparseLU<SparseMatrix<double>, COLAMDOrdering<int>>(A);

    /** forJay
     * We convert A * x = b to the equivalent system B * y = c, where B = P^T * A * P, y = P^T * x, c = P^T * b
     * In other words, the new RHS for the linear solver is c = P^T * b
     * The linear solver gives y. Left-multiply by P to get x (i.e., x = P * y)
    **/

    x = perm.transpose() * b;
    lu.matrixLU().template triangularView<UnitLower>().solveInPlace(x);
    lu.matrixLU().template triangularView<Upper>().solveInPlace(x);

#endif

    // Convert back to chemical_state
    ChemicalState result = {};
    for (int i = 0; i < n_variables; ++i)
    {
        result[i] = x(i);
    }

    return result;
}


ChemicalState
direct_solve(const SparseMatrix<double>& A, const ChemicalState& b_)
{
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;

    return
    direct_solve(A, b_, solver, true);
}

#else

template<typename T>
void swap_gen(T& a, T& b)
{
    T temp = a;
    a = b;
    b = temp;
}

SpeciesJacobian invert_jacobian(const SpeciesJacobian& J)
{
    SpeciesJacobian A = J; // make a copy
    SpeciesJacobian inv = {};

    // Initialize inv to the identity matrix
    for (int i = 0; i < n_variables; ++i)
        inv[i][i] = 1.0;

    for (int i = 0; i < n_variables; ++i)
    {
        // Pivot: find the max row in column i
        int max_row = i;
        for (int k = i + 1; k < n_variables; ++k)
        {
            if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
                max_row = k;
        }

        // Swap rows in both A and inv
        swap_gen(A[i], A[max_row]);
        swap_gen(inv[i], inv[max_row]);

        double pivot = A[i][i];
        if (std::abs(pivot) < 1e-14)
        {
            std::cerr << "Matrix is singular or nearly singular.\n";
            return SpeciesJacobian{}; // or throw exception
        }

        // Normalize the pivot row
        for (int j = 0; j < n_variables; ++j)
        {
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        // Eliminate column i in other rows
        for (int k = 0; k < n_variables; ++k)
        {
            if (k == i) continue;
            double factor = A[k][i];
            for (int j = 0; j < n_variables; ++j)
            {
                A[k][j] -= factor * A[i][j];
                inv[k][j] -= factor * inv[i][j];
            }
        }
    }

    return inv;
}

ChemicalState direct_solve(const SpeciesJacobian& J, const ChemicalState& rhs)
{
    return invert_jacobian(J) * rhs;
}
#endif