#ifdef CHEMGEN_EIGEN


ChemicalState
direct_solve(const SparseMatrix<double>& A, const ChemicalState& b_, SparseLU<SparseMatrix<double>, NaturalOrdering<int>>& lu, bool compute)
{
    Matrix<double, n_variables, 1> b(n_variables);

    for (int i = 0; i < n_variables; ++i)
    {
        b(i) = b_[i];
    }

    Matrix<double, n_variables, 1> x;

    SparseMatrix<double> B;

    // Permutations
    PermutationMatrix<n_variables, n_variables> perm; // column perm (input)
    perm.indices() = perm_indices;
    B = perm.transpose() * A * perm;

    PermutationMatrix<n_variables, n_variables> row_perm; // note: this corresponds to PA = LU

#ifdef CHEMGEN_EIGEN_SPARSE
    MatrixXd A_, B_;

    lu.setPivotThreshold(0.);

    if (compute) lu.compute(B);
    if (lu.info() != Success) std::cerr << "lu factorization failed!" << std::endl;
    x = perm * lu.solve(perm.transpose() * b);

    row_perm.indices() = lu.rowsPermutation().indices();

#if 0
    A_ = A; B_ = B;
    std::cout << "\n\nA = \n" << A_ << std::endl;
    std::cout << "\n\nB = \n" << B_ << std::endl;
    std::cout << "\n\nlu.rowsPermutation().indices() = \n" << lu.rowsPermutation().indices() << std::endl;
    // std::exit(0);
#endif

#else

    if (compute) lu = SparseLU<SparseMatrix<double>, NaturalOrdering<int>>(B);

    /** forJay
     * We convert A * x = b to the equivalent system B * y = c, where B = P^T * A * P, y = P^T * x, c = P^T * b
     * In other words, the new RHS for the linear solver is c = P^T * b
     * The linear solver gives y. Left-multiply by P to get x (i.e., x = P * y)
    **/

    x = perm.transpose() * b;
    lu.matrixLU().template triangularView<UnitLower>().solveInPlace(x);
    lu.matrixLU().template triangularView<Upper>().solveInPlace(x);
    x = perm * x;

    row_perm.indices() = lu.permutationP().indices();

    // Check expected sparsity - can't easily do this for SparseLU since can't directly access L and U matrices
    // Need extremely small dt to prevent pivoting since can't directly control pivoting threshold (unlike for Sparse LU)
    MatrixXd m_lu = lu.matrixLU(); // LU decomposition where L (excluding unit diagonal) is stored in lower triangular part and U is stored in upper triangular part

#if 1
    for (int i = 0; i < n_variables; ++i)
    {
        for (int j = 0; j < n_variables; ++j)
        {
            if (m_lu(i,j) != 0. && sp_lu_perm[i][j] == 0)
            {
                std::cout << "i = " << i << std::endl;
                std::cout << "j = " << j << std::endl;
                std::cout << "m_lu(i,j) = " << m_lu(i,j) << std::endl;

                std::exit(0);
            }
        }
    }
#endif

#endif

#if 1
    // Check that row permutation is identity
    for (int i = 0; i < n_variables; ++i)
    {
        if (row_perm.indices()(i) != i)
        {
            std::cout << "i = " << i << std::endl;
            std::cout << "row_perm.indices()(i) = " << row_perm.indices()(i) << std::endl;

            std::exit(0);
        }
    }
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
    SparseLU<SparseMatrix<double>, NaturalOrdering<int>> solver;

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