ChemicalState 
backwards_euler(ChemicalState y,  
                const double& dt,
                //...........................................
                std::chrono::duration<double>& NN_total_time,
                std::chrono::duration<double>& P_total_time,
                int& cvs_iter,
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

        int n = 0;
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

            // #define CHEMGEN_DIRECT_SOLVER
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dy = invert_jacobian(A) * res;
            //.............................
            // NOT WORKING (DO NOT USE)
            // # define CHEMGEN_DIRECT_PINN
            #elif CHEMGEN_DIRECT_PINN
            Species dy = train_pinn(A,res);
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            #else
            Species dy = gmres_solve(A, res, 
                                    NN_total_time,
                                    P_total_time,
                                    cvs_iter,
                                    tol);
            #endif

            y_guess = y_guess + dy;
            n = iter;
            if (norm2(dy) < 1e-10) // check convergence
            {
                // return set_chemical_state(y[0], y_guess); 
                break;
            };
        }

        // std::cout << "# Newton iterations = " << n + 1 << std::endl;
        
        return set_chemical_state(y[0], y_guess); //sets energy to zero to signal integration broke
}