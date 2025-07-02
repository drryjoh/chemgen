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