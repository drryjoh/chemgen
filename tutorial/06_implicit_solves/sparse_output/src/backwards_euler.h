
ChemicalState 
backwards_euler(ChemicalState y,  
                const double& dt,
                // int & csv_iter, 
                int  &n_gmres_iter_total,
                double abs_newton_tol = newton_abs_tol_default,
                double rel_newton_tol = newton_rel_tol_default,
                int  max_iter = newton_iter_max_default)
{
        ChemicalState y_init = y;
        ChemicalState y_guess = y;

        SparseMatrix<double> A;
        GMRES<SparseMatrix<double>, Preconditioner> solver;
        
        

        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            double temperature_ = temperature(y[0], y_guess);
            ChemicalState f = source(y_guess, temperature_);
            ChemicalState res = {};

            
            {
                A = source_jacobian(y_guess, temperature_, -1, 1/dt);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(A);
            // std::ostringstream ss;
            // ss << "BE_DATA/jacobian_" << csv_iter << ".csv"; //change path
            // std::ofstream out(ss.str());
            // for (int i = 0; i < denseA.rows(); ++i) 
            // {
            //     for (int j = 0; j < denseA.cols(); ++j) 
            //     {
            //         out << denseA(i, j);
            //         if (j + 1 < denseA.cols())
            //             out << ",";
            //     }
            //     out << "\n";
            // }
            // csv_iter += 1;
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            //RHS
            res = scale_gen(-double(1)/dt, y_guess - y_init) + f;
            //Solve
            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dy = direct_solve(A, res, solver, true);
            #else
            ChemicalState dy = gmres_solve(A, res, solver, true, n_gmres_iter_total);
            #endif

            //Increment
            y_guess = y_guess + dy;

            if (error_norm(dy, y_guess, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                // std::cout << "# Newton iterations = " << iter + 1 << std::endl;

#if defined(CHEMGEN_EIGEN) && !defined(CHEMGEN_EIGEN_SPARSE)
                if (max_iter > 10)
                {
                    int n_nonzeros = 0;

                    for (int i = 0; i < n_variables; ++i)
                    {
                        for (int j = 0; j < n_variables; ++j)
                        {
                            if (A(i,j) != 0) n_nonzeros++;
                        }
                    }

                    // std::cout << "A = \n" << A << std::endl;
                    std::cout << "n_nonzeros = \n" << n_nonzeros << std::endl;
                }
#endif
                return y_guess;
            }
        }

        // std::cout << "Newton solve did not converge!" << std::endl;

        return y_guess; //sets energy to zero to signal integration broke
}