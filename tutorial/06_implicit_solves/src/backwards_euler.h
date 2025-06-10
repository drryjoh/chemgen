ChemicalState 
backwards_euler(ChemicalState y,  
                const double& dt,
                //...........................................
                // bool final_step,
                // int cvs_iter,
                // std::chrono::duration<double>& NN_total_time,
                // std::chrono::duration<double>& P_total_time,
                //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                double tol = 1e-12, 
                int max_iter = 10
                ) 
{        
/*******************************************************************************
1.  -   initialize 
            y^n_k
    -   use NEWTON-RAPHSON METHOD to get temperature guess

    -   initialize species guess 
            y^{n+1}_k

2.  -   get temperature at current guess 
            T @ y^{n+1}_k

    -   get source term at current species guess and temperature 
            S(y^{n+1}_k)

    -   get jacobian of source term at current species guess and temperature  
            dS(y^{n+1}_k) / y^{n+1}_k

    -   intialize a identity matrix for the jacobian

3.  LHS:
    -   we want to solve the residual such that 
            R(y^{n+1}_k) = 0

    -   thus we taylor expand the residual around 
            y^{n_k} 
        and set it to zero

    -   R(y^{n+1}_k) = (dR/dy^{n+1})@k * (y^{n+1}_k - y^{n_k}) 
                     = 0

    -   thus we can rewrite the jacobian of the residual as 
            (dR/dy^{n+1})@k = I/dt - dS(y^{n+1}_k) / y^{n+1}_k

    -   thus our LHS is 
            [I/dt - (dS(y^{n+1}_k) / y^{n+1}_k)] * (y^{n+1}_k - y^{n_k})

4. RHS:
    -   our RHS based of the LHS is then 
            R(y^{n+1}_k) 
        which is defined as 
            - 1/dt * (y^{n+1}_k - y^{n_k}) + S(y^{n+1}_k) 
        from above

5. Solve:
    -   our total equation is now 
            (dR/dy^{n+1})@k * (y^{n+1}_k - y^{n_k}) = R(y^{n+1}_k)
        or alternatively 
            (dR/dy^{n+1})@k * dy = R(y^{n+1}_k)
        or to be exact 
            [I/dt - (dS(y^{n+1}_k) / y^{n+1}_k)] * dy = - 1/dt * (y^{n+1}_k - y^{n_k}) + S(y^{n+1}_k)

    -   we use a linear solver (GMRES) to approx dy

6.  -   update y using dy
            y^{n+1}_{k+1} = y^{n+1}_k + dy

    -   check for convergence
***********************************************************************************/

        Species y_init = get_species(y);
        double temperature_guess = temperature(y[0], y_init);
        Species y_guess = get_species(y);

        //...............................................
        // ////////////////////////////////////////////
        // //// THIS IS USED FOR TRAINING PURPOSES ////
        // ////////////////////////////////////////////
        // SpeciesJacobian last_A;
        // Species last_res;
        // Species last_dy;
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        for (int iter = 0; iter < 5; ++iter) 
        {
            double temperature_ = temperature(y[0], y_guess); 
            Species f = source(y_guess, temperature_);
            SpeciesJacobian J = source_jacobian(y_guess, temperature_); 
            SpeciesJacobian A = jacobian_I();
            Species res = {};
            
            for (int i = 0; i < n_species; ++i) A[i][i] = A[i][i]/dt; 
            A = A - J; 

            res = scale_gen(-double(1)/dt, y_guess - y_init) + f; 


            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dy = invert_jacobian(A) * res;
            //.............................
            // NOT WORKING (DO NOT USE)
            #elif CHEMGEN_DIRECT_PINN
            Species dy = train_pinn(A,res);
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            #else
            Species dy = gmres_solve(A, res, tol);
            #endif

            //.............................................
            // ////////////////////////////////////////////
            // //// THIS IS USED FOR TRAINING PURPOSES ////
            // ////////////////////////////////////////////
            // last_A = invert_jacobian(A);
            // last_res = res;
            // last_dy = dy;
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            y_guess = y_guess + dy;

            if (norm2(dy) < 1e-10) // check convergence
            {
                // return set_chemical_state(y[0], y_guess); 
                break;
            };
        }

        //............................................................................
        // ////////////////////////////////////////////
        // //// THIS IS USED FOR TRAINING PURPOSES ////
        // ////////////////////////////////////////////
        // if (final_step)
        // {
        //     std::string output_dir = "./neural_net/data_all/"; 

        //     // std::ofstream A_file(output_dir + "A_inv_{cvs_iter}.csv");
        //     std::ofstream A_file(output_dir + "A_inv_" + std::to_string(cvs_iter) + ".csv");
        //     if (A_file.is_open())
        //     {
        //         for (int i = 0; i < n_species; ++i)
        //         {
        //             for (int j = 0; j < n_species; ++j)
        //             {
        //                 A_file << last_A[i][j];
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
            
        //     // std::ofstream res_file(output_dir + "res.csv");
        //     std::ofstream res_file(output_dir + "res_" + std::to_string(cvs_iter) + ".csv");
        //     if (res_file.is_open())
        //     {
        //         for (int j = 0; j < n_species; ++j)
        //         {
        //             res_file << last_res[j];
        //             if (j != n_species - 1)
        //                 res_file << ",";
        //         }
        //         res_file << "\n";
        //         res_file.close();
        //     }
        //     else
        //     {
        //         std::cerr << "Error opening file for residual output!" << std::endl;
        //     }

        //     // std::ofstream dy_file(output_dir + "dy.csv");
        //     std::ofstream dy_file(output_dir + "dy_" + std::to_string(cvs_iter) + ".csv");
        //     if (dy_file.is_open())
        //     {
        //         for (int j = 0; j < n_species; ++j)
        //         {
        //             dy_file << last_dy[j];
        //             if (j != n_species - 1)
        //                 dy_file << ",";
        //         }
        //         dy_file << "\n";
        //         dy_file.close();
        //     }
        //     else
        //     {
        //         std::cerr << "Error opening file for dy output!" << std::endl;
        //     }
        // }
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        return set_chemical_state(y[0], y_guess); //sets energy to zero to signal integration broke
}