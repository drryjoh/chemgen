ChemicalState 
backwards_euler(ChemicalState y,  
                const double& dt,
                double tol, 
                int max_iter,
                //=====================================================================
                ////////////////////////////////////////////
                //// THIS IS USED FOR TRAINING PURPOSES ////
                ////////////////////////////////////////////
                // bool final_step,
                std::chrono::duration<double>& NN_total_time,
                std::chrono::duration<double>& P_total_time
                //=====================================================================
                ) 
{        
        // initialize y^n_k
        Species y_init = get_species(y);

        // NEWTON-RAPHSON METHOD to get temperature guess
        double temperature_guess = temperature(y[0], y_init);

        // initialize species guess y^{n+1}_k
        Species y_guess = get_species(y);

        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        SpeciesJacobian J_init = source_jacobian(y_guess, temperature_guess); 
        SpeciesJacobian A_init = jacobian_I();
        for (int i = 0; i < n_species; ++i) A_init[i][i] = A_init[i][i]/dt; 
            A_init = A_init - J_init; 

        auto NN_start = std::chrono::high_resolution_clock::now();

        SpeciesJacobian P = cnn_2(A_init);

        auto NN_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> NN_duration = NN_end - NN_start;
        // std::cout << "[NN Part] Time elapsed: " << NN_duration.count() << " seconds" << std::endl;
        NN_total_time += NN_duration;
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        //=====================================================================
        // ////////////////////////////////////////////
        // //// THIS IS USED FOR TRAINING PURPOSES ////
        // ////////////////////////////////////////////
        // SpeciesJacobian last_A;
        // Species last_res;
        // Species last_dy;
        //===================================================================

        /* 
        - solve dy/dt - S(y) = 0 
        - using 5 newton iterations (not time-stepping)
        - resiudal = - 1/dt * (y^{n+1}_k - y^{n_k}) + S(y^{n+1}_k) = R
        */
        for (int iter = 0; iter < 5; ++iter) 
        {
            // get temperature at current guess y^{n+1}_k
            double temperature_ = temperature(y[0], y_guess); 

            // get source term at current species guess and temperature S(y^{n+1}_k)
            Species f = source(y_guess, temperature_);

            // get jacobian of source term at current species guess and temperature  
            // dS(y^{n+1}_k) / y^{n+1}_k
            SpeciesJacobian J = source_jacobian(y_guess, temperature_); 

            // intialize a identity matrix for the jacobian
            SpeciesJacobian A = jacobian_I();

            // allocate/define variable for residual
            Species res = {};
            
            //LHS
            /*
            - we want to solve the residual such that R(y^{n+1}_k) = 0
            - thus we taylor expand the residual around y^{n_k} and set it to zero
            - R(y^{n+1}_k) = (dR/dy^{n+1})@k * (y^{n+1}_k - y^{n_k}) = 0
            - thus we can rewrite the jacobian of the residual as (dR/dy^{n+1})@k = I/dt - dS(y^{n+1}_k) / y^{n+1}_k
            - thus our LHS is [I/dt - (dS(y^{n+1}_k) / y^{n+1}_k)] * (y^{n+1}_k - y^{n_k})
            */
            for (int i = 0; i < n_species; ++i) A[i][i] = A[i][i]/dt; 
            A = A - J; 

            //RHS
            // our RHS based of the LHS is then R(y^{n+1}_k) which is defined as: - 1/dt * (y^{n+1}_k - y^{n_k}) + S(y^{n+1}_k) from above
            res = scale_gen(-double(1)/dt, y_guess - y_init) + f; 

            //Solve
            /* 
            - our total equation is now: (dR/dy^{n+1})@k * (y^{n+1}_k - y^{n_k}) = R(y^{n+1}_k)
            - or alternatively: (dR/dy^{n+1})@k * dy = R(y^{n+1}_k)
            - or to be exact: [I/dt - (dS(y^{n+1}_k) / y^{n+1}_k)] * dy = - 1/dt * (y^{n+1}_k - y^{n_k}) + S(y^{n+1}_k)
            - we use a linear solver (GMRES) to approx dy
            */
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            auto P_start = std::chrono::high_resolution_clock::now();

            A = operator*(P, A);
            res = operator*(P, res);

            auto P_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> P_duration = P_end - P_start;
            // std::cout << "[Preconditioning Part] Time elapsed: " << P_duration.count() << " seconds" << std::endl;
            P_total_time += P_duration;
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dy = invert_jacobian(A) * res;
            #else
            Species dy = gmres_solve(A, res, tol = 1e-10);
            #endif

            //=====================================================================
            // ////////////////////////////////////////////
            // //// THIS IS USED FOR TRAINING PURPOSES ////
            // ////////////////////////////////////////////
            // last_A = invert_jacobian(A);
            // last_res = res;
            // last_dy = dy;
            //===================================================================

            //Increment
            // y^{n+1}_{k+1} = y^{n+1}_k + dy
            y_guess = y_guess + dy;

            // check convergence
            if (norm2(dy) < 1e-10)
            {

                //=====================================================================
                // ////////////////////////////////////////////
                // //// THIS IS USED FOR TRAINING PURPOSES ////
                // ////////////////////////////////////////////
                // if (final_step)
                // {
                //     std::string output_dir = "./neural_net/data/"; 

                //     // save A
                //     std::ofstream A_file(output_dir + "A_inv.csv");
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
                    
                //     //save r
                //     std::ofstream res_file(output_dir + "res.csv");
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

                //     // Save dy (1D vector)
                //     std::ofstream dy_file(output_dir + "dy.csv");
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
                //=====================================================================

                return set_chemical_state(y[0], y_guess); 
            };
        }
        
        return set_chemical_state(y[0], y_guess); //sets energy to zero to signal integration broke
}