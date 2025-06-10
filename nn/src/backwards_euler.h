
ChemicalState 
backwards_euler(ChemicalState y,  
                const double& dt, 
                double tol = 1e-10, 
                int max_iter = 10) 
{
        
        
        Species y_init = get_species(y);
        double temperature_guess = temperature(y[0], y_init);
        Species y_guess = get_species(y);
        for (int iter = 0; iter < 5; ++iter) 
        {
            double temperature_ = temperature(y[0], y_guess);
            Species f = source(y_guess, temperature_);
            SpeciesJacobian J = source_jacobian(y_guess, temperature_);
            SpeciesJacobian A = jacobian_I();
            Species res = {};
            
            //LHS
            for (int i = 0; i < n_species; ++i) A[i][i] = A[i][i]/dt; 
            A = A - J;

            //RHS
            res = scale_gen(-double(1)/dt, y_guess - y_init) + f;
            //Solve
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dy = invert_jacobian(A) * res;
            #else
            Species dy = gmres_solve(A, res, tol = 1e-10);
            #endif

            //Increment
            y_guess = y_guess + dy;

            if (norm2(dy) < 1e-10)
            {
                return set_chemical_state(y[0], y_guess); 
            };
        }
        return set_chemical_state(y[0], y_guess); //sets energy to zero to signal integration broke
}