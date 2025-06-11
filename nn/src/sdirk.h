/*
SpeciesJacobian M_inv = inverse_diagonal(J);
Species preconditioned_res = apply_diagonal(M_inv, res);
SpeciesJacobian preconditioned_J = apply_diagonal(M_inv, J);
*/


ChemicalState 
sdirk2(ChemicalState y,  
       const double& dt, 
       double tol = 1e-10,
       int max_iter = 10) 
{
        Species y_init = get_species(y);
        double gamma = 1.0 - 1.0 / sqrt_gen(double(2));
        double one_minus_gamma = 1.0 - gamma;
        double temperature_init = temperature(y[0], y_init);
        Species k1 = source(y_init, temperature_init);
        Species k2 = k1;
        const SpeciesJacobian I = jacobian_I();

        //stage1
        for (int iter = 0; iter < 5; ++iter) 
        {
            Species y_stage = y_init + scale_gen(gamma * dt, k1);
            double temperature_ = temperature(y[0], y_stage);
            Species f_val = source(y_stage, temperature_);
            Species res = f_val - k1;
            if (norm2(res) < 1e-6)
            {
                break;
            };
            SpeciesJacobian J = I -scale_gen(gamma * dt, source_jacobian(y_stage, temperature_));
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dk = invert_jacobian(J) * res;
            #else
            Species dk = gmres_solve(J, res);
            #endif
            k1 = k1 + dk;
        }

        //stage2
        for (int iter = 0; iter < 5; ++iter) 
        {
            //y_stage = yn + dt * (one_minus_gamma * k1 + gamma * k2)
            Species y_stage = y_init + scale_gen(one_minus_gamma * dt, k1) + scale_gen(gamma * dt, k2);
            double temperature_ = temperature(y[0], y_stage);
            Species f_val = source(y_stage, temperature_);
            Species res = f_val - k2;
            if (norm2(res) < 1e-6)
            {
                break;
            };
            SpeciesJacobian J = I -scale_gen(gamma * dt, source_jacobian(y_stage, temperature_));
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dk = invert_jacobian(J) * res;
            #else
            Species dk = gmres_solve(J, res);
            #endif
            k2 = k2 + dk;
        }

        return 
        set_chemical_state(y[0], 
                           y_init + scale_gen(dt * one_minus_gamma, k1) + scale_gen(gamma * dt, k2));
}



ChemicalState 
sdirk4(ChemicalState y,  
       const double& dt, 
       double tol = 1e-10,
       int max_iter = 10) 
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

        Species y_init = get_species(y);
        double temperature_init = temperature(y[0], y_init);
        Species k1 = source(y_init, temperature_init);
        Species k2 = k1;
        Species k3 = k1;
        Species k4 = k1;
        Species k5 = k1;
        const SpeciesJacobian I = jacobian_I();
/*
        Confirmed Stage 1
*/

        for (int iter = 0; iter < 5; ++iter) 
        {
            Species y_stage = y_init + scale_gen(gamma * dt, k1);
            double temperature_ = temperature(y[0], y_stage);
            Species f_val = source(y_stage, temperature_);
            Species res = f_val - k1;
            if (norm2(res) < 1e-6)
            {
                break;
            };
            SpeciesJacobian J = I -scale_gen(gamma * dt, source_jacobian(y_stage, temperature_));
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dk = invert_jacobian(J) * res;
            #else
            Species dk = gmres_solve(J, res);
            #endif
            k1 = k1 + dk;
        }

        /*
        # Stage 2 confirmed
        */
        k2 = k1;
        for (int iter = 0; iter < 5; ++iter) 
        {
            //y_stage = yn + dt * (one_minus_gamma * k1 + gamma * k2)
            Species y_stage = y_init + scale_gen(dt, scale_gen(a21, k1) + scale_gen(gamma, k2));
            double temperature_ = temperature(y[0], y_stage);
            Species f_val = source(y_stage, temperature_);
            Species res = f_val - k2;
            if (norm2(res) < 1e-6)
            {
                break;
            };
            SpeciesJacobian J = I -scale_gen(gamma * dt, source_jacobian(y_stage, temperature_));
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dk = invert_jacobian(J) * res;
            #else
            Species dk = gmres_solve(J, res);
            #endif
            k2 = k2 + dk;
        }
        /*
        #stage 3 
        */
        k3 = k2;
        for (int iter = 0; iter < 5; ++iter) 
        {
            //y_stage = yn + dt * (a31 * k[0] + a32 * k[1] + gamma * k[2])
            Species y_stage = y_init + scale_gen(dt, scale_gen(a31, k1) + scale_gen(a32, k2) + scale_gen(gamma, k3));
            double temperature_ = temperature(y[0], y_stage);
            Species f_val = source(y_stage, temperature_);
            Species res = f_val - k3;
            if (norm2(res) < 1e-6)
            {
                break;
            };
            SpeciesJacobian J = I -scale_gen(gamma * dt, source_jacobian(y_stage, temperature_));
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dk = invert_jacobian(J) * res;
            #else
            Species dk = gmres_solve(J, res);
            #endif
            k3 = k3 + dk;
        }

        /*
        stage 4
        */
        k4 = k3;
        for (int iter = 0; iter < 5; ++iter) 
        {
            //y_stage = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2] + gamma * k[3])
            Species y_stage = y_init +
            scale_gen(dt, scale_gen(a41, k1) + scale_gen(a42, k2) + scale_gen(a43, k3) + scale_gen(gamma, k4));
            double temperature_ = temperature(y[0], y_stage);
            Species f_val = source(y_stage, temperature_);
            Species res = f_val - k4;
            if (norm2(res) < 1e-6)
            {
                break;
            };

            SpeciesJacobian J = I - scale_gen(gamma * dt, source_jacobian(y_stage, temperature_));
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dk = invert_jacobian(J) * res;
            #else
            Species dk = gmres_solve(J, res);
            #endif
            k4 = k4 + dk;
        }

        /*
        stage 5
        */
        k5 = k4;
        for (int iter = 0; iter < 5; ++iter) 
        {
            //y_stage = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2] + gamma * k[3])
            Species y_stage = y_init + 
            scale_gen(dt, scale_gen(a51, k1) + scale_gen(a52, k2) + scale_gen(a53, k3) + scale_gen(a54, k4)+ scale_gen(gamma, k5));
            double temperature_ = temperature(y[0], y_stage);
            Species f_val = source(y_stage, temperature_);
            Species res = f_val - k5;
            if (norm2(res) < 1e-6)
            {
                break;
            };
            SpeciesJacobian J = I -scale_gen(gamma * dt, source_jacobian(y_stage, temperature_));
            #ifdef CHEMGEN_DIRECT_SOLVER
            Species dk = invert_jacobian(J) * res;
            #else
            Species dk = gmres_solve(J, res);
            #endif


            k5 = k5 + dk;
        }


        return 
        set_chemical_state(y[0], 
                           y_init + scale_gen(dt, (scale_gen(b1, k1) + scale_gen(b2, k2) + scale_gen(b3, k3) + scale_gen(b4, k4) + scale_gen(b5, k5))));
}