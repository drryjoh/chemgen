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

// #define CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL
#ifdef CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL

Species apply_gauss_seidel(const SpeciesJacobian &A, const Species &v)
{
    Species z = {};
    for (int i = 0; i < n_species; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
        {
            sum += A[i][j] * z[j]; // Forward substitution
        }

        double diag = A[i][i];
        z[i] = (v[i] - sum) / (std::abs(diag) > 1e-14 ? diag : 1.0);
    }
    return z;
}
#endif

// ..................................................................
// #define CHEMGEN_PRECONDITIONER_NN
#ifdef CHEMGEN_PRECONDITIONER_NN

std::array<double, 96 * 96> flatten_matrix(SpeciesJacobian A)
{
    std::array<double, 96 * 96> flat_A;
    for (int i = 0; i < 96; i++)
    {
        for (int j = 0; j < 96; j++)
        {
            flat_A[(i*96)+j] = A[i][j];
        }
    }
    return flat_A;
}

SpeciesJacobian unflatten_matrix(std::array<double, 96 * 96> A)
{
    SpeciesJacobian Q;
    for (int i = 0; i < 96; i++)
    {
        for (int j = 0; j < 96; j++)
        {
            Q[i][j] = A[(i*96)+j];
        }
    }
    return Q;
}

std::pair<std::array<double, 96 * 96>, std::array<double, 96 * 96>> split_matrix(std::array<double, 96 * 96 * 2> LU)
{
    std::array<double, 96 * 96> L;
    std::array<double, 96 * 96> U;
    for (int i = 0; i < 96; i++)
    {
        for (int j = 0; j < 96; j++)
        {
            L[(i*96)+j] = LU[(i*96*2)+j];
            U[(i*96)+j] = LU[(i*96*2)+j+96];
        }
    }
    return {L,U};
}
#endif

// #define CHEMGEN_PRECONDITIONER_NN_DIAGONAL
#ifdef CHEMGEN_PRECONDITIONER_NN_DIAGONAL

std::array<double, 96 * 96> flatten_matrix(SpeciesJacobian A)
{
    std::array<double, 96 * 96> flat_A;
    for (int i = 0; i < 96; i++)
    {
        for (int j = 0; j < 96; j++)
        {
            flat_A[(i*96)+j] = A[i][j];
        }
    }
    return flat_A;
}

SpeciesJacobian unflatten_matrix(std::array<double, 96 * 96> A)
{
    SpeciesJacobian Q;
    for (int i = 0; i < 96; i++)
    {
        for (int j = 0; j < 96; j++)
        {
            Q[i][j] = A[(i*96)+j];
        }
    }
    return Q;
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
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

















//////////////////
//// ORIGINAL ////
//////////////////
Species gmres_solve(const SpeciesJacobian &A, const Species &b,
                    double abs_tol = 1e-8, double rel_tol = 1e-4)
{
    int max_iter = 100;
    Species x = {};
    Species cs = {};
    Species sn = {};
    bool GS = false;
#if defined(CHEMGEN_PRECONDITIONER_JACOBI)
    SpeciesJacobian P = inverse_diagonal(A);
    SpeciesJacobian A_ = apply_diagonal(P, A);
    Species b_ = apply_diagonal(P, b);
    Species r = b_ - (A_ * x);
#elif defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
    SpeciesJacobian A_ = A;
    Species b_ = apply_gauss_seidel(A, b);
    GS = true;
    Species r = b_ - (A_ * x);
#else
    SpeciesJacobian A_ = A;
    Species b_ = b;
    Species r = b_ - (A_ * x);
#endif

    double norm2_r = norm2(r);
    double norm2_A = norm2(A_);

    if (norm2_r < abs_tol)
    {
        return x;
    }

    std::array<Species, n_species + 1> V = {};
    V[0] = scale_gen(inv_gen(norm2_r), r);

    std::array<std::array<double, n_species + 1>, n_species + 1> H = {};
    std::array<double, n_species + 1> g = {};
    g[0] = norm2(r);

    int final_iter = 0;

    for (int j = 0; j < n_species; ++j)
    {
        final_iter = j;
        
        // Apply gauss-seidel if yes
        Species w;
        #if defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
        {
            w = apply_gauss_seidel(A, A * V[j]);
        }
        #else
        {
            w = A_ * V[j];
        }
        #endif

        // Modified Gram-Schmidt
        for (int i = 0; i <= j; ++i)
        {
            H[i][j] = dot(V[i], w);
            w = w - H[i][j] * V[i];
        }

        H[j + 1][j] = norm2(w);
        if (H[j + 1][j] < abs_tol * norm2_A)
            break;
        V[j + 1] = scale_gen(inv_gen(H[j + 1][j]), w);

        // Apply Givens rotations to new column of H
        for (int i = 0; i < j; ++i)
        {
            double temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
            H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
            H[i][j] = temp;
        }
        // Compute new Givens rotation
        double a = H[j][j];
        double b_h = H[j + 1][j];
        double r_val = std::sqrt(a * a + b_h * b_h);
        cs[j] = a / r_val;
        sn[j] = b_h / r_val;

        // Apply to H and g
        H[j][j] = r_val;
        H[j + 1][j] = 0.0;

        g[j + 1] = 0.0; // Ensure valid memory before rotation
        double temp_g = cs[j] * g[j] + sn[j] * g[j + 1];
        g[j + 1] = -sn[j] * g[j] + cs[j] * g[j + 1];
        g[j] = temp_g;

        // Convergence check
        double res_norm = std::abs(g[j + 1]);
        if (res_norm < abs_tol || res_norm < rel_tol * norm2_r)
            break;
    }

    // Solve least squares problem Hy = g using back-substitution on H (upper
    // triangular approx)
    std::array<double, n_species> y = {};
    for (int i = final_iter; i >= 0; --i)
    {
        double sum = 0;
        for (int j = i + 1; j <= final_iter; ++j)
            sum += H[i][j] * y[j];
        y[i] = (g[i] - sum) / H[i][i];
    }
    Species result = {};
    for (int i = 0; i < n_species; ++i)
        result = result + y[i] * V[i];

    return result;
}
















/////////////////
//// TIMING /////
/////////////////
Species gmres_solve(const SpeciesJacobian &A, const Species &b, 
                    std::chrono::duration<double>& NN_total_time,
                    std::chrono::duration<double>& P_total_time,
                    int& cvs_iter,
                    double abs_tol = 1e-8, double rel_tol = 1e-4)
{
    int max_iter = 100;
    Species x = {};
    Species cs = {};
    Species sn = {};
#if defined(CHEMGEN_PRECONDITIONER_JACOBI)
    SpeciesJacobian P = inverse_diagonal(A);
    // Species r = P*(b - A*x);
    Species r = apply_diagonal(P, b - A*x);
    SpeciesJacobian A_ = A;
#elif defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
    SpeciesJacobian A_ = A;
    Species b_ = apply_gauss_seidel(A, b);
    GS = true;
    Species r = b_ - (A_ * x);
//.............................................
#elif defined(CHEMGEN_PRECONDITIONER_NN)
    ///////////////////////
    //// REGULAR A_inv ////
    ///////////////////////
    auto p_start = std::chrono::high_resolution_clock::now(); //start preconditioning timer
    std::array<double, 96 * 96> A_flat = flatten_matrix(A);
    auto nn_start = std::chrono::high_resolution_clock::now(); //start inference timer
    std::array<double, 96 * 96> P_flat;
    if (cvs_iter < 340) 
    {
        P_flat = MLP_2_BE(A_flat);
    } 
    else 
    {
        P_flat = MLP_2_BE_second(A_flat);
    }
    // SpeciesJacobian P = CNN_1(A); //CNN
    auto nn_end = std::chrono::high_resolution_clock::now(); //end inference timer
    SpeciesJacobian P = unflatten_matrix(P_flat);
    Species r = P*(b - A*x);
    SpeciesJacobian A_ = A;
    auto p_end = std::chrono::high_resolution_clock::now(); //end preconditioning timer
    std::chrono::duration<double> nn_duration = nn_end - nn_start;
    NN_total_time += nn_duration; //total inference time
    std::chrono::duration<double> p_duration = p_end - p_start - nn_duration;
    P_total_time += p_duration; //total precon time

    //////////////////////////
    //// LU DECOMPOSITION ////
    //////////////////////////
    // auto p_start = std::chrono::high_resolution_clock::now(); //start preconditioning timer
    // std::array<double, 96 * 96> A_flat = flatten_matrix(A);
    // auto nn_start = std::chrono::high_resolution_clock::now(); //start inference timer
    // std::array<double, 96 * 96 * 2> LU_flat = MLP_x_1(A_flat);
    // auto nn_end = std::chrono::high_resolution_clock::now(); //end inference timer
    // std::chrono::duration<double> nn_duration = nn_end - nn_start;
    // NN_total_time += nn_duration; //total inference time
    // auto [L_flat, U_flat] = split_matrix(LU_flat);
    // SpeciesJacobian L = unflatten_matrix(L_flat);
    // SpeciesJacobian U = unflatten_matrix(U_flat);
    // // SpeciesJacobian P = operator*(L, U);
    // Species r = L * ( U * (b - A * x));
    // SpeciesJacobian A_ = A;
    // auto p_end = std::chrono::high_resolution_clock::now(); //end preconditioning timer
    // std::chrono::duration<double> p_duration = p_end - p_start - nn_duration;
    // P_total_time += p_duration; //total precon time
#elif defined(CHEMGEN_PRECONDITIONER_NN_DIAGONAL)
    std::array<double, 96 * 96> A_flat = flatten_matrix(A);
    std::array<double, 96 * 96> P_flat = TINIER(A_flat);
    SpeciesJacobian P = unflatten_matrix(P_flat);
    Species r = apply_diagonal(P, b - A*x);
    SpeciesJacobian A_ = A;
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#else
    SpeciesJacobian A_ = A;
    Species b_ = b;
    Species r = b_ - (A_ * x);
#endif

    double norm2_r = norm2(r);
    // double norm2_A = norm2(A_);
    // double norm2_A = norm2(apply_diagonal(P, A_));

    if (norm2_r < abs_tol)
    {
        // std::cout << "# GMRES iterations = 0" << std::endl;
        return x;
    }

    std::array<Species, n_species + 1> V = {};
    V[0] = scale_gen(inv_gen(norm2_r), r);

    std::array<std::array<double, n_species + 1>, n_species + 1> H = {};
    std::array<double, n_species + 1> g = {};
    g[0] = norm2(r);

    int final_iter = 0;

    for (int j = 0; j < n_species; ++j)
    {
        final_iter = j;
        
        // Apply gauss-seidel if yes
        Species w;
        #if defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
        {
            w = apply_gauss_seidel(A, A * V[j]);
        }
        #elif defined(CHEMGEN_PRECONDITIONER_NN)
        {
            w = P * (A * V[j]);
            // w = L * (U * (A * V[j]));
        }
        #elif defined(CHEMGEN_PRECONDITIONER_JACOBI)
        {
            // w = P * (A * V[j]);
            w = apply_diagonal(P, A * V[j]);
        }
        #else
        {
            w = A_ * V[j];
        }
        #endif

        // Modified Gram-Schmidt
        for (int i = 0; i <= j; ++i)
        {
            H[i][j] = dot(V[i], w);
            w = w - H[i][j] * V[i];
        }

        H[j + 1][j] = norm2(w);
        // if (H[j + 1][j] < abs_tol * norm2_A)
        //     break;
        V[j + 1] = scale_gen(inv_gen(H[j + 1][j]), w);

        // Apply Givens rotations to new column of H
        for (int i = 0; i < j; ++i)
        {
            double temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
            H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
            H[i][j] = temp;
        }
        // Compute new Givens rotation
        double a = H[j][j];
        double b_h = H[j + 1][j];
        double r_val = std::sqrt(a * a + b_h * b_h);
        cs[j] = a / r_val;
        sn[j] = b_h / r_val;

        // Apply to H and g
        H[j][j] = r_val;
        H[j + 1][j] = 0.0;

        g[j + 1] = 0.0; // Ensure valid memory before rotation
        double temp_g = cs[j] * g[j] + sn[j] * g[j + 1];
        g[j + 1] = -sn[j] * g[j] + cs[j] * g[j + 1];
        g[j] = temp_g;

        // Convergence check
        double res_norm = std::abs(g[j + 1]);
        if (res_norm < abs_tol || res_norm < rel_tol * norm2_r)
            break;
    }

    // Solve least squares problem Hy = g using back-substitution on H (upper
    // triangular approx)
    std::array<double, n_species> y = {};
    for (int i = final_iter; i >= 0; --i)
    {
        double sum = 0;
        for (int j = i + 1; j <= final_iter; ++j)
            sum += H[i][j] * y[j];
        y[i] = (g[i] - sum) / H[i][i];
    }
    Species result = {};
    for (int i = 0; i < n_species; ++i)
        result = result + y[i] * V[i];

    // std::cout << "# GMRES iterations = " << final_iter + 1 << std::endl;
    // std::cout << "ITER #: " << cvs_iter << std::endl;
    // if (cvs_iter == 340){
    //     std::cout << "========================" << std::endl;
    //     std::cout << "========================" << std::endl;
    //     std::cout << "========================" << std::endl;
    // }
    // if (cvs_iter == 400){
    //     std::cout << "================" << std::endl;
    // }
    cvs_iter += 1;

    return result;
}
























// /////////////////////
// //// SAVING DATA ////
// /////////////////////
// Species gmres_solve(const SpeciesJacobian &A, const Species &b,
//                     int &cvs_iter,
//                     double abs_tol = 1e-8, double rel_tol = 1e-4)
// {

//     int max_iter = 100;
//     Species x = {};
//     Species cs = {};
//     Species sn = {};

// #if defined(CHEMGEN_PRECONDITIONER_JACOBI)
//     SpeciesJacobian P = inverse_diagonal(A);
//     // Species r = P*(b - A*x);
//     Species r = apply_diagonal(P, b - A*x);
//     SpeciesJacobian A_ = A;
// #elif defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
//     SpeciesJacobian A_ = A;
//     Species b_ = apply_gauss_seidel(A, b);
//     GS = true;
//     Species r = b_ - (A_ * x);
// //.............................................
// #elif defined(CHEMGEN_PRECONDITIONER_NN)
//     //// REGULAR A_inv ////
//     // std::array<double, 96 * 96> A_flat = flatten_matrix(A);
//     // std::array<double, 96 * 96> P_flat = MLP_2(A_flat);
//     // // SpeciesJacobian P = CNN_1(A);
//     // SpeciesJacobian P = unflatten_matrix(P_flat);
//     // Species r = P*(b - A*x);
//     // SpeciesJacobian A_ = A;

//     //// LU DECOMPOSITION ////
//     std::array<double, 96 * 96> A_flat = flatten_matrix(A);
//     std::array<double, 96 * 96 * 2> LU_flat = MLP_LU_1(A_flat);
//     auto [L_flat, U_flat] = split_matrix(LU_flat);
//     SpeciesJacobian L = unflatten_matrix(L_flat);
//     SpeciesJacobian U = unflatten_matrix(U_flat);
//     // SpeciesJacobian P = operator*(L, U);
//     Species r = L * ( U * (b - A * x));
//     SpeciesJacobian A_ = A;
// #elif defined(CHEMGEN_PRECONDITIONER_NN_DIAGONAL)
//     std::array<double, 96 * 96> A_flat = flatten_matrix(A);
//     std::array<double, 96 * 96> P_flat = TINIER(A_flat);
//     SpeciesJacobian P = unflatten_matrix(P_flat);
//     Species r = apply_diagonal(P, b - A*x);
//     SpeciesJacobian A_ = A;
// //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// #else
//     SpeciesJacobian A_ = A;
//     Species b_ = b;
//     Species r = b_ - (A_ * x);
// #endif

//     double norm2_r = norm2(r);

//     if (norm2_r < abs_tol)
//     {
//         // std::cout << "# GMRES iterations = 0" << std::endl;
//         return x;
//     }

//     std::array<Species, n_species + 1> V = {};
//     V[0] = scale_gen(inv_gen(norm2_r), r);

//     std::array<std::array<double, n_species + 1>, n_species + 1> H = {};
//     std::array<double, n_species + 1> g = {};
//     g[0] = norm2(r);

//     int final_iter = 0;

//     for (int j = 0; j < n_species; ++j)
//     {
//         final_iter = j;
        
//         // Apply gauss-seidel if yes
//         Species w;
//         #if defined(CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL)
//         {
//             w = apply_gauss_seidel(A, A * V[j]);
//         }
//         #elif defined(CHEMGEN_PRECONDITIONER_NN)
//         {
//             // w = P * (A * V[j]);
//             w = L * (U * (A * V[j]));
//         }
//         #elif defined(CHEMGEN_PRECONDITIONER_JACOBI)
//         {
//             // w = P * (A * V[j]);
//             w = apply_diagonal(P, A * V[j]);
//         }
//         #else
//         {
//             w = A_ * V[j];
//         }
//         #endif

//         // Modified Gram-Schmidt
//         for (int i = 0; i <= j; ++i)
//         {
//             H[i][j] = dot(V[i], w);
//             w = w - H[i][j] * V[i];
//         }

//         H[j + 1][j] = norm2(w);
//         // if (H[j + 1][j] < abs_tol * norm2_A)
//         //     break;
//         V[j + 1] = scale_gen(inv_gen(H[j + 1][j]), w);

//         // Apply Givens rotations to new column of H
//         for (int i = 0; i < j; ++i)
//         {
//             double temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
//             H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
//             H[i][j] = temp;
//         }
//         // Compute new Givens rotation
//         double a = H[j][j];
//         double b_h = H[j + 1][j];
//         double r_val = std::sqrt(a * a + b_h * b_h);
//         cs[j] = a / r_val;
//         sn[j] = b_h / r_val;

//         // Apply to H and g
//         H[j][j] = r_val;
//         H[j + 1][j] = 0.0;

//         g[j + 1] = 0.0; // Ensure valid memory before rotation
//         double temp_g = cs[j] * g[j] + sn[j] * g[j + 1];
//         g[j + 1] = -sn[j] * g[j] + cs[j] * g[j + 1];
//         g[j] = temp_g;

//         // Convergence check
//         double res_norm = std::abs(g[j + 1]);
//         if (res_norm < abs_tol || res_norm < rel_tol * norm2_r)
//             break;
//     }

//     // Solve least squares problem Hy = g using back-substitution on H (upper
//     // triangular approx)
//     std::array<double, n_species> y = {};
//     for (int i = final_iter; i >= 0; --i)
//     {
//         double sum = 0;
//         for (int j = i + 1; j <= final_iter; ++j)
//             sum += H[i][j] * y[j];
//         y[i] = (g[i] - sum) / H[i][i];
//     }
//     Species result = {};
//     for (int i = 0; i < n_species; ++i)
//         result = result + y[i] * V[i];

//     if (final_iter + 1 < n_species)
//     {
//         // std::cout << "# GMRES iterations = " << final_iter + 1 << std::endl;
//     }

//     //............................................................................
//     ////////////////////////////////////////////
//     //// THIS IS USED FOR TRAINING PURPOSES ////
//     ////////////////////////////////////////////
//     std::string output_dir = "./neural_net/SDIRK2_DATA/"; 

//     // std::ofstream A_file(output_dir + "A_inv_{cvs_iter}.csv");
//     std::ofstream A_file(output_dir + "A_" + std::to_string(cvs_iter) + ".csv");
//     if (A_file.is_open())
//     {
//         for (int i = 0; i < n_species; ++i)
//         {
//             for (int j = 0; j < n_species; ++j)
//             {
//                 A_file << A[i][j];
//                 if (j != n_species - 1)
//                     A_file << ","; // builds reach row
//             }
//             A_file << "\n";
//         }
//         A_file.close();
//     }
//     else
//     {
//         std::cerr << "Error opening file for A output!" << std::endl;
//     }

//     SpeciesJacobian invA = invert_jacobian(A);
//     std::ofstream invA_file(output_dir + "A_inv_" + std::to_string(cvs_iter) + ".csv");
//     if (invA_file.is_open())
//     {
//         for (int i = 0; i < n_species; ++i)
//         {
//             for (int j = 0; j < n_species; ++j)
//             {
//                 invA_file << invA[i][j];
//                 if (j != n_species - 1)
//                     invA_file << ","; // builds reach row
//             }
//             invA_file << "\n";
//         }
//         invA_file.close();
//     }
//     else
//     {
//         std::cerr << "Error opening file for A output!" << std::endl;
//     }
    
//     std::ofstream b_file(output_dir + "b_" + std::to_string(cvs_iter) + ".csv");
//     if (b_file.is_open())
//     {
//         for (int j = 0; j < n_species; ++j)
//         {
//             b_file << b[j];
//             if (j != n_species - 1)
//                 b_file << ",";
//         }
//         b_file << "\n";
//         b_file.close();
//     }
//     else
//     {
//         std::cerr << "Error opening file for residual output!" << std::endl;
//     }

//     std::ofstream x_file(output_dir + "x_" + std::to_string(cvs_iter) + ".csv");
//     if (x_file.is_open())
//     {
//         for (int j = 0; j < n_species; ++j)
//         {
//             x_file << result[j];
//             if (j != n_species - 1)
//                 x_file << ",";
//         }
//         x_file << "\n";
//         x_file.close();
//     }
//     else
//     {
//         std::cerr << "Error opening file for dy output!" << std::endl;
//     }
//     //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//     cvs_iter += 1;

//     return result;
// }