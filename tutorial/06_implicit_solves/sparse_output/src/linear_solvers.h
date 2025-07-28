#define CHEMGEN_PRECONDITIONER_NN

#define CHEMGEN_EIGEN

#ifdef CHEMGEN_PRECONDITIONER_JACOBI

ChemicalState inverse_diagonal(SpeciesJacobian J) 
{
    ChemicalState diag_inv;
    for (int i = 0; i < n_variables; ++i)
    {
        diag_inv[i] = std::abs(J[i][i]) > 1e-14 ? 1.0 / J[i][i] : 0.0;
    }
    return diag_inv;
}
#endif

#ifdef CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL

ChemicalState
apply_gauss_seidel(SpeciesJacobian A, ChemicalState v)
{
    ChemicalState z = {};
    for (int i = 0; i < n_variables; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
            sum += A[i][j] * z[j]; // Already computed values of z

        // Diagonal + lower
        double diag = A[i][i];
        z[i] = (v[i] - sum) / (std::abs(diag) > 100.*std::numeric_limits<double>::epsilon() ? diag : 1.0); // Avoid div by zero
    }
    return z;
}
#endif

#ifdef CHEMGEN_PRECONDITIONER_NN

std::array<double, (n_species + 1)*(n_species + 1)> flatten_jacobian(const SpeciesJacobian& A)
{
    std::array<double, (n_species + 1)*(n_species + 1)> flat_A;
    for (int i = 0; i < (n_species + 1); i++)
    {
        for (int j = 0; j < (n_species + 1); i++)
        {
            flat_A[(i*(n_species + 1))+j] = A[i][j];
        }
    }
    return flat_A;
}


const SpeciesJacobian& unflatten_jacobian(const std::array<double, (n_species + 1)*(n_species + 1)>& flat_A)
{
    SpeciesJacobian Q;
    for (int i = 0; i < (n_species + 1); ++i) {
        for (int j = 0; j < (n_species + 1); ++j) {
            Q[i][j] = flat_A[i*(n_species + 1) + j];
        }
    }
    return Q;
}
#endif

#ifdef CHEMGEN_EIGEN
#include "./nn_preconditioner.hpp"
#include "./custom_preconditioners_eigen.h"

using Preconditioner =

#ifdef CHEMGEN_PRECONDITIONER_JACOBI
DiagonalPreconditioner<double>;

#elif defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
GaussSeidelPreconditioner<double, SparseMatrix<double>>;

#elif defined(CHEMGEN_PRECONDITIONER_ILU)
IncompleteLUT<double>;

#elif defined(CHEMGEN_PRECONDITIONER_NN)
NNPreconditioner<double>;

#else
IdentityPreconditioner;

#endif


ChemicalState
gmres_solve(const SparseMatrix<double>& A, const ChemicalState& b_, GMRES<SparseMatrix<double>, Preconditioner>& gmres, bool compute, int &n_gmres_iter_total, double abs_tol=linear_abs_tol_default, double rel_tol=linear_rel_tol_default)
{
    Matrix<double, n_variables, 1> b(n_variables);

    for (int i = 0; i < n_variables; ++i)
    {
        b(i) = b_[i];
    }
    gmres.setMaxIterations(n_variables); // total number of iterations including restarts
    gmres.set_restart(n_variables); // # iterations before restart (Full-size Krylov subspace)

    if (compute)
    {
#if defined(CHEMGEN_PRECONDITIONER_ILU)

        // gmres.preconditioner().setDroptol(NumTraits<double>::dummy_precision()); // default value
        gmres.preconditioner().setDroptol(0.0001); //0.0001 is most optimal

        //gmres.preconditioner().setFillfactor(10); // default value
        gmres.preconditioner().setFillfactor(n_variables * n_variables); // allow for maximum fill-in, such that only droptol controls ilu sparsity

#endif

        gmres.compute(A); // preconditioner initialized here

        if (gmres.preconditioner().info() != Success) std::cerr << "Preconditioner failed!" << std::endl;
    }

    // eigen uses relative tolerance - https://eigen.tuxfamily.org/dox/unsupported/GMRES_8h_source.html
    double tol = std::max(rel_tol, divide(abs_tol, gmres.preconditioner().solve(b).norm()));
    gmres.setTolerance(tol);

    Matrix<double, n_variables, 1> x = gmres.solve(b);

    // std::cout << "gmres.iterations() = " << gmres.iterations() << std::endl;
    if (gmres.info() != Success) std::cerr << "GMRES failed!" << std::endl;

    n_gmres_iter_total += gmres.iterations();

    // Convert back to chemical_state
    ChemicalState result = {};
    for (int i = 0; i < n_variables; ++i)
    {
        result[i] = x(i);
    }

    return result;
}


ChemicalState
gmres_solve(const SparseMatrix<double>& A, const ChemicalState& b_, int &n_gmres_iter_total, double abs_tol=linear_abs_tol_default, double rel_tol=linear_rel_tol_default)
{
    GMRES<SparseMatrix<double>, Preconditioner> solver;

    return
    gmres_solve(A, b_, solver, true, n_gmres_iter_total, abs_tol, rel_tol);
}

#else


ChemicalState
gmres_solve(const SpeciesJacobian& A, const ChemicalState& b, int &n_gmres_iter_total, double abs_tol=linear_abs_tol_default, double rel_tol=linear_rel_tol_default)
{
    int max_iter = 100;
    ChemicalState x = {};
    ChemicalState cs = {};
    ChemicalState sn = {};

    // Note: if initial guess for x is nonzero, need to replace b below with (b - A*x)
#ifdef CHEMGEN_PRECONDITIONER_JACOBI
    ChemicalState diag_inv = inverse_diagonal(A);
    ChemicalState r = diag_inv * (b);
#elif defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
    ChemicalState r = apply_gauss_seidel(A, b);
#elif CHEMGEN_PRECONDITIONER_NN
    SpeciesJacobian A_ = A;
    std::array<double, n_variables * n_variables> A_flat = flatten_jacobian(A);
    std::array<double, n_variables * n_variables> P_flat = MLP_BE(A_flat);
    SpeciesJacobian P = unflatten_jacobian(P_flat);
    Species r = P*(b - A_*x); 
#else
    ChemicalState r = b;
#endif

    double norm2_r = norm2(r);

    if (norm2_r < std::numeric_limits<double>::epsilon())
    {
        // std::cout << "# GMRES iterations = 0" << std::endl;
        return x;
    }

    std::array<ChemicalState, n_variables + 1> V = {};
    V[0] = scale_gen(inv_gen(norm2_r), r);

    std::array<std::array<double, n_variables>, n_variables + 1> H = {};
    std::array<double, n_variables + 1> g = {};
    g[0] = norm2(r);

    int final_iter = 0;

    for (int j = 0; j < n_variables; ++j)
    {
        final_iter = j;

        ChemicalState w = A * V[j];
#ifdef CHEMGEN_PRECONDITIONER_JACOBI
        w = diag_inv * w;
#elif defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
        w = apply_gauss_seidel(A, w);
#elif defined(CHEMGEN_PRECONDITIONER_NN)
        w = P * w;
#endif

        #ifdef CHEMGEN_PRECONDITIONER_NN
        {
            w = P * (A * V[j]);
        }
        #endif
        // Modified Gram-Schmidt
        for (int i = 0; i <= j; ++i)
        {
            H[i][j] = dot(V[i], w);
            w = w - H[i][j] * V[i];
        }

        H[j+1][j] = norm2(w);

        V[j+1] = scale_gen(inv_gen(H[j+1][j]),w); // not used if j == n_variables - 1

        // Apply Givens rotations to new column of H
        for (int i = 0; i < j; ++i)
        {
            double temp = cs[i] * H[i][j] + sn[i] * H[i+1][j];
            H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
            H[i][j] = temp;
        }
        // Compute new Givens rotation
        double a = H[j][j];
        double b_h = H[j+1][j];
        double r_val = std::sqrt(a * a + b_h * b_h);
        cs[j] = a / r_val;
        sn[j] = b_h / r_val;

        // Apply to H and g
        H[j][j] = r_val;
        H[j+1][j] = 0.0;
        
        double temp_g = cs[j] * g[j];
        g[j+1] = -sn[j] * g[j];
        g[j] = temp_g;
        
        // Convergence check
        double res_norm = std::abs(g[j+1]);
        if (res_norm < abs_tol || res_norm < rel_tol * norm2_r)
            break;
    }

    // Solve least squares problem Hy = g using back-substitution on H (upper triangular, rectangular)
    std::array<double, n_variables> y = {};
    for (int i = final_iter; i >= 0; --i)
    {
        double sum = 0;
        for (int j = i + 1; j <= final_iter; ++j)
            sum += H[i][j] * y[j];
        y[i] = (g[i] - sum) / H[i][i];
    }

    ChemicalState result = {};
    for (int i = 0; i < final_iter + 1; ++i)
        result = result + y[i] * V[i];

    if (final_iter + 1 < n_variables)
    {
        // std::cout << "# GMRES iterations = " << final_iter + 1 << std::endl;
    }

    n_gmres_iter_total += final_iter + 1;

    return result;
}
#endif