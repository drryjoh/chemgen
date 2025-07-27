/*
SpeciesJacobian M_inv = inverse_diagonal(J);
ChemicalState preconditioned_res = apply_diagonal(M_inv, res);
SpeciesJacobian preconditioned_J = apply_diagonal(M_inv, J);
*/


ChemicalState
sdirk2(ChemicalState y,
       const double& dt,
       //int & csv_iter,
       int &n_gmres_iter_total,
       double abs_newton_tol = newton_abs_tol_default,
       double rel_newton_tol = newton_rel_tol_default,
       int  max_iter = newton_iter_max_default)
{
        ChemicalState y_init = y;
        double gamma = 1.0 - 1.0 / sqrt_gen(double(2));
        double one_minus_gamma = 1.0 - gamma;
        ChemicalState k1 = source(y_init);
        ChemicalState k2 = k1;

        SparseMatrix<double> J;
        GMRES<SparseMatrix<double>, Preconditioner> solver;
        
        

        //stage1
        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            ChemicalState y_stage = y_init + scale_gen(gamma * dt, k1);
            double temperature_ = temperature(y[0], y_stage);

            ChemicalState f_val = source(y_stage, temperature_);
            ChemicalState res = f_val - k1;

            
            {
                J = source_jacobian(y_stage, temperature_, -gamma*dt, 1);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(J);
            // std::ostringstream ss;
            // ss << "SD2_DATA/jacobian_" << csv_iter << ".csv"; //change path
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

            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dk = direct_solve(J, res, solver, true);
            #else
            ChemicalState dk = gmres_solve(J, res, solver, true, n_gmres_iter_total);
            #endif
            k1 = k1 + dk;
            if (error_norm(dk, k1, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                break;
            };
        }
        

        //stage2
        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            //y_stage = yn + dt * (one_minus_gamma * k1 + gamma * k2)
            ChemicalState y_stage = y_init + scale_gen(one_minus_gamma * dt, k1) + scale_gen(gamma * dt, k2);
            double temperature_ = temperature(y[0], y_stage);
            ChemicalState f_val = source(y_stage, temperature_);
            ChemicalState res = f_val - k2;

            
            {
                J = source_jacobian(y_stage, temperature_, -gamma*dt, 1);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(J);
            // std::ostringstream ss;
            // ss << "SD2_DATA/jacobian_" << csv_iter << ".csv"; //change path
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

            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dk = direct_solve(J, res, solver, true);
            #else
            ChemicalState dk = gmres_solve(J, res, solver, true, n_gmres_iter_total);
            #endif
            k2 = k2 + dk;
            if (error_norm(dk, k2, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                break;
            };
        }

        return
        y_init + scale_gen(dt * one_minus_gamma, k1) + scale_gen(gamma * dt, k2);
}



ChemicalState
sdirk4(ChemicalState y,
       const double& dt,
       //int & csv_iter,
       int &n_gmres_iter_total,
       double abs_newton_tol = newton_abs_tol_default,
       double rel_newton_tol = newton_rel_tol_default,
       int  max_iter = newton_iter_max_default)
{
        double gamma = .26666666666666666666666666666666670;
        double a11 = .26666666666666666666666666666666670;
        double a21 = .50000000000000000000000000000000000;
        double a22 = .26666666666666666666666666666666670;
        double a31 = .35415395284327323162274618585298200;
        double a32 = -.5415395284327323162274618585298197e-1;
        double a33 = .26666666666666666666666666666666670;
        double a41 = .8515494131138652076337791881433756e-1;
        double a42 = -.6484332287891555171683963466229754e-1;
        double a43 = .7915325296404206392428857585141242e-1;
        double a44 = .26666666666666666666666666666666670;
        double a51 = 2.1001157005669327779706120559990740;
        double a52 = -.76778002844459768133431021850622760;
        double a53 = 2.3998163610800263980947462052738800;
        double a54 = -2.9988186998690281613977147094333940;
        double a55 = .26666666666666666666666666666666670;

        double b1   = 2.1001157005669327779706120559990740;
        double b2   = -.76778002844459768133431021850622760;
        double b3   = 2.3998163610800263980947462052738800;
        double b4   = -2.9988186998690281613977147094333940;
        double b5   = .26666666666666666666666666666666670;

        ChemicalState y_init = y;
        ChemicalState k1 = source(y_init);
        ChemicalState k2 = k1;
        ChemicalState k3 = k1;
        ChemicalState k4 = k1;
        ChemicalState k5 = k1;

        SparseMatrix<double> J;
        GMRES<SparseMatrix<double>, Preconditioner> solver;
        
        
/*
        Confirmed Stage 1
*/

        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            ChemicalState y_stage = y_init + scale_gen(gamma * dt, k1);
            double temperature_ = temperature(y[0], y_stage);
            ChemicalState f_val = source(y_stage, temperature_);
            ChemicalState res = f_val - k1;

            
            {
                J = source_jacobian(y_stage, temperature_, -gamma*dt, 1);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(J);
            // std::ostringstream ss;
            // ss << "SD4_DATA/jacobian_" << csv_iter << ".csv"; //change path
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

            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dk = direct_solve(J, res, solver, true);
            #else
            ChemicalState dk =  gmres_solve(J, res, solver, true, n_gmres_iter_total);
            #endif
            k1 = k1 + dk;
            if (error_norm(dk, k1, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                break;
            };
        }
        

        /*
        # Stage 2 confirmed
        */
        k2 = k1;
        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            //y_stage = yn + dt * (one_minus_gamma * k1 + gamma * k2)
            ChemicalState y_stage = y_init + scale_gen(dt, scale_gen(a21, k1) + scale_gen(gamma, k2));
            double temperature_ = temperature(y[0], y_stage);
            ChemicalState f_val = source(y_stage, temperature_);
            ChemicalState res = f_val - k2;

            
            {
                J = source_jacobian(y_stage, temperature_, -gamma*dt, 1);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(J);
            // std::ostringstream ss;
            // ss << "SD4_DATA/jacobian_" << csv_iter << ".csv"; //change path
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

            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dk = direct_solve(J, res, solver, true);
            #else
            ChemicalState dk = gmres_solve(J, res, solver, true, n_gmres_iter_total);
            #endif
            k2 = k2 + dk;
            if (error_norm(dk, k2, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                break;
            };
        }
        
        /*
        #stage 3
        */
        k3 = k2;
        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            //y_stage = yn + dt * (a31 * k[0] + a32 * k[1] + gamma * k[2])
            ChemicalState y_stage = y_init + scale_gen(dt, scale_gen(a31, k1) + scale_gen(a32, k2) + scale_gen(gamma, k3));
            double temperature_ = temperature(y[0], y_stage);
            ChemicalState f_val = source(y_stage, temperature_);
            ChemicalState res = f_val - k3;

            
            {
                J = source_jacobian(y_stage, temperature_, -gamma*dt, 1);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(J);
            // std::ostringstream ss;
            // ss << "SD4_DATA/jacobian_" << csv_iter << ".csv"; //change path
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

            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dk = direct_solve(J, res, solver, true);
            #else
            ChemicalState dk = gmres_solve(J, res, solver, true, n_gmres_iter_total);
            #endif
            k3 = k3 + dk;
            if (error_norm(dk, k3, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                break;
            };
        }
        

        /*
        stage 4
        */
        k4 = k3;
        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            //y_stage = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2] + gamma * k[3])
            ChemicalState y_stage = y_init +
            scale_gen(dt, scale_gen(a41, k1) + scale_gen(a42, k2) + scale_gen(a43, k3) + scale_gen(gamma, k4));
            double temperature_ = temperature(y[0], y_stage);
            ChemicalState f_val = source(y_stage, temperature_);
            ChemicalState res = f_val - k4;

            
            {
                J = source_jacobian(y_stage, temperature_, -gamma*dt, 1);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(J);
            // std::ostringstream ss;
            // ss << "SD4_DATA/jacobian_" << csv_iter << ".csv"; //change path
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

            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dk = direct_solve(J, res, solver, true);
            #else
            ChemicalState dk = gmres_solve(J, res, solver, true, n_gmres_iter_total);
            #endif
            k4 = k4 + dk;
            if (error_norm(dk, k4, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                break;
            };
        }
        

        /*
        stage 5
        */
        k5 = k4;
        for (int iter = 0; iter < max_iter; ++iter)
        {
            
            //y_stage = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2] + gamma * k[3])
            ChemicalState y_stage = y_init +
            scale_gen(dt, scale_gen(a51, k1) + scale_gen(a52, k2) + scale_gen(a53, k3) + scale_gen(a54, k4)+ scale_gen(gamma, k5));
            double temperature_ = temperature(y[0], y_stage);
            ChemicalState f_val = source(y_stage, temperature_);
            ChemicalState res = f_val - k5;;

            
            {
                J = source_jacobian(y_stage, temperature_, -gamma*dt, 1);
            }
            //................................................................
            // Eigen::MatrixXd denseA = Eigen::MatrixXd(J);
            // std::ostringstream ss;
            // ss << "SD4_DATA/jacobian_" << csv_iter << ".csv"; //change path
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

            #ifdef CHEMGEN_DIRECT_SOLVER
            ChemicalState dk = direct_solve(J, res, solver, true);
            #else
            ChemicalState dk = gmres_solve(J, res, solver, true, n_gmres_iter_total);
            #endif


            k5 = k5 + dk;
            if (error_norm(dk, k5, abs_newton_tol, rel_newton_tol) < 1.0)
            {
                break;
            };
        }

        return
        y_init + scale_gen(dt, (scale_gen(b1, k1) + scale_gen(b2, k2) + scale_gen(b3, k3) + scale_gen(b4, k4) + scale_gen(b5, k5)));
}