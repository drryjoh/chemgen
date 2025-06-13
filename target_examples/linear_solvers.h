#define CHEMGEN_PRECONDITIONER_JACOBI


#ifdef CHEMGEN_PRECONDITIONER_JACOBI

SpeciesJacobian inverse_diagonal(SpeciesJacobian J) 
{
    SpeciesJacobian M_inv = {};
    for (int i = 0; i < n_species; ++i)
    {
        M_inv[i][i] = std::abs(J[i][i]) > 1e-14 ? 1.0 / J[i][i] : 0.0;
    }
    return M_inv;
}


SpeciesJacobian apply_diagonal(SpeciesJacobian P, SpeciesJacobian A) 
{
    SpeciesJacobian M = {};
    for (int i = 0; i < n_species; ++i)
    {
        for (int j = 0; j < n_species; ++j)
        {
            M[i][j] = P[i][i] * A[i][j];
        }
    }
    return M;
}


Species apply_diagonal(SpeciesJacobian P, Species b) 
{
    Species m = {};
    for (int i = 0; i < n_species; ++i)
    {
        m[i] = P[i][i] * b[i];
    }
    return m;
}
#endif

#ifdef CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL

Species
apply_gauss_seidel(SpeciesJacobian A, Species v)
{
    Species z = {};
    for (int i = 0; i < n_species; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
            sum += A[i][j] * z[j]; // Already computed values of z

        // Diagonal + lower
        double diag = A[i][i];
        z[i] = (v[i] - sum) / (std::abs(diag) > 1e-14 ? diag : 1.0); // Avoid div by zero
    }
    return z;
}
#endif


Species
gmres_solve(const SpeciesJacobian& A, const Species& b, double abs_tol = 1e-8, double rel_tol = 1e-4) 
{
    int max_iter = 100;
    Species x = {};
    Species cs = {};
    Species sn = {};
#ifdef CHEMGEN_PRECONDITIONER_JACOBI
    SpeciesJacobian P = inverse_diagonal(A);
    SpeciesJacobian A_ = apply_diagonal(P, A);
    Species b_ = apply_diagonal(P, b);
#else
    SpeciesJacobian A_ = A;
    Species b_ = b;
#endif

    Species r = b_ - (A_ * x);
    double norm2_r = norm2(r);
    double norm2_A = norm2(A_);
    
    if (norm2_r < abs_tol)
    {
        return x;
    }

    std::array<Species, n_species + 1> V = {};
    V[0] = scale_gen(inv_gen(norm2_r), r); // stores orthonormal basis vectors of Krylov subspace

    std::array<std::array<double, n_species>, n_species + 1> H = {}; // Hessenberg matrix that's modified in place at each iteration to yield the upper triangular (rectangular) matrix
    std::array<double, n_species + 1> g = {};
    g[0] = norm2(r);

    int final_iter = 0;

    for (int j = 0; j < n_species; ++j)
    {
        final_iter = j;
        Species w = A_ * V[j];

        // Modified Gram-Schmidt
        for (int i = 0; i <= j; ++i)
        {
            H[i][j] = dot(V[i], w);
            w = w - H[i][j] * V[i];
        }

        H[j+1][j] = norm2(w);
        if (H[j+1][j] < abs_tol * norm2_A)
            break;
        V[j+1] = scale_gen(inv_gen(H[j+1][j]),w); // last one not used

        // Apply all previous Givens rotations to new column of H
        // This "new column of H" is the new column of the Hessenberg matrix
        // Each Givens rotation matrix is orthogonal; product of orthogonal matrices is also orthogonal, which is "Q" in the QR decomposition of the Hessenberg matrix
        for (int i = 0; i < j; ++i)
        {
            double temp = cs[i] * H[i][j] + sn[i] * H[i+1][j];
            H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
            H[i][j] = temp;
        }
        // At this point, H is almost a triangular matrix - only non-triangular portion is the single nonzero first-subdiagonal entry in column j

        // Compute new Givens rotation
        // Note that the Givens rotation matrix defined on Wikipedia (https://en.wikipedia.org/wiki/Generalized_minimal_residual_method#Solving_the_least_squares_problem)
        // is misleading since the rotation is applied to the new column j, whereas Wikipedia says it's always applied to the nth (last) column
        // In other words, when it says H_n is an (n+1)-by-n matrix, it should actually be something like (n+1)-by-(j+1),
        // and the cs and sn entries of the rotation matrix should be placed accordingly
        double a = H[j][j];
        double b_h = H[j+1][j];
        double r_val = std::sqrt(a * a + b_h * b_h);
        cs[j] = a / r_val;
        sn[j] = b_h / r_val;

        // Apply to H and g
        // H is modified in place to yield the upper triangular matrix
        // This can be easily done since successive Hessenberg matrices differ only by a row of zeros and a column
        H[j][j] = r_val;
        H[j+1][j] = 0.0;
        
        // Apply Givens rotation to g (one at a time)
        // Note that applying the orthogonal matrix of the QR decomposition of the Hessenberg matrix is equivalent to applying all Givens rotations successively
        g[j+1] = 0.0; // Ensure valid memory before rotation
        double temp_g = cs[j] * g[j] + sn[j] * g[j+1]; // g[j+1] is zero
        g[j+1] = -sn[j] * g[j] + cs[j] * g[j+1]; // g[j+1] is zero
        g[j] = temp_g;
        
        // Convergence check
        double res_norm = std::abs(g[j+1]); // this is the component of g that cannot be zeroed out by H*y since the row (j+1) of H is all zero
        if (res_norm < abs_tol || res_norm < rel_tol * norm2_r)
            break;
    }


    // Solve least squares problem Hy = g using back-substitution on H (upper triangular, rectangular)
    std::array<double, n_species> y = {};
    for (int i = final_iter; i >= 0; --i)
    {
        double sum = 0;
        for (int j = i + 1; j <= final_iter; ++j)
            sum += H[i][j] * y[j];
        y[i] = (g[i] - sum) / H[i][i];

    }
    Species result = {};
    for (int i = 0; i < n_species; ++i) // can optimize this
        result = result + y[i] * V[i];

    return result;
}