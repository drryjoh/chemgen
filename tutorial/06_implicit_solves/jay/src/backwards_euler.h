
ChemicalState 
backwards_euler(ChemicalState y,  
                const double& dt,
                int  &n_gmres_iter_total,
                double abs_newton_tol = newton_abs_tol_default,
                double rel_newton_tol = newton_rel_tol_default,
                int  max_iter = newton_iter_max_default)
{
        ChemicalState y_init = y;
        ChemicalState y_guess = y;
        for (int iter = 0; iter < max_iter; ++iter)
        {
            double temperature_ = temperature(y[0], y_guess);
            ChemicalState f = source(y_guess, temperature_);
            SparseMatrix<double> A = source_jacobian(y_guess, temperature_, -1, 1/dt);
            ChemicalState res = {};

            //RHS
            res = scale_gen(-double(1)/dt, y_guess - y_init) + f;
            //Solve
            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dy = invert_jacobian(A) * res;
            #else
            ChemicalState dy = gmres_solve(A, res, n_gmres_iter_total);
            #endif

            //Increment
            y_guess = y_guess + dy;

            if (error_norm(dy, y_guess, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                // std::cout << "# Newton iterations = " << iter + 1 << std::endl;
                return y_guess;
            }
        }

        // std::cout << "Newton solve did not converge!" << std::endl;

        return y_guess; //sets energy to zero to signal integration broke
}