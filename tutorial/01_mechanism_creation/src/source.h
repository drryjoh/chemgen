
        
        Species source(const Species& species, const double& temperature)  
        {
            Species net_production_rates = {double(0)};
            double inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature);
            double log_temperature = log_gen(temperature);
            Reactions gibbs_reactions = gibbs_reaction(log_temperature);
            double pressure_ = pressure(species, temperature);
            double mixture_concentration = pressure_ * inv_universal_gas_constant_temperature;
        double forward_reaction_0 = call_forward_reaction_0(temperature, log_temperature);
        double forward_reaction_1 = call_forward_reaction_1(temperature, log_temperature);
        double forward_reaction_2 =  call_forward_reaction_2(species, temperature, log_temperature, mixture_concentration);
        double forward_reaction_3 =  call_forward_reaction_3(temperature, pressure_);
        double forward_reaction_4 =  call_forward_reaction_4(species, temperature, log_temperature, mixture_concentration);
        double forward_reaction_5 =  call_forward_reaction_5(species, temperature, log_temperature, mixture_concentration);
        double forward_reaction_6 =  call_forward_reaction_6(species, temperature, log_temperature, mixture_concentration);
        double equilibrium_constant_0 = multiply(exp_gen(-gibbs_reactions[0]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
        double rate_of_progress_0 = multiply(pow_gen2(species[1]) * species[3], forward_reaction_0) - multiply(pow_gen2(species[5]), divide(forward_reaction_0, equilibrium_constant_0));
        double rate_of_progress_1 = multiply(pow_gen2(species[1]) * species[3], forward_reaction_1);
        double equilibrium_constant_2 = multiply(exp_gen(-gibbs_reactions[2]), pow_gen2(pressure_atmosphere() * inv_universal_gas_constant_temperature));
        double rate_of_progress_2 = multiply(species[5], forward_reaction_2) - multiply(pow_gen2(species[0]) * species[2], divide(forward_reaction_2, equilibrium_constant_2));
        double equilibrium_constant_3 = multiply(exp_gen(-gibbs_reactions[3]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
        double rate_of_progress_3 = multiply(pow_gen2(species[1]) * species[3], forward_reaction_3) - multiply(pow_gen2(species[5]), divide(forward_reaction_3, equilibrium_constant_3));
        double equilibrium_constant_4 = multiply(exp_gen(-gibbs_reactions[4]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
        double rate_of_progress_4 = multiply(pow_gen2(species[1]) * species[3], forward_reaction_4) - multiply(pow_gen2(species[5]), divide(forward_reaction_4, equilibrium_constant_4));
        double equilibrium_constant_5 = multiply(exp_gen(-gibbs_reactions[5]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
        double rate_of_progress_5 = multiply(pow_gen2(species[1]) * species[3], forward_reaction_5) - multiply(pow_gen2(species[5]), divide(forward_reaction_5, equilibrium_constant_5));
        double equilibrium_constant_6 = multiply(exp_gen(-gibbs_reactions[6]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
        double rate_of_progress_6 = multiply(pow_gen2(species[1]) * species[3], forward_reaction_6) - multiply(pow_gen2(species[5]), divide(forward_reaction_6, equilibrium_constant_6));

        net_production_rates[0] = double(2.0) * rate_of_progress_2;
        net_production_rates[1] = double(-2.0) * rate_of_progress_0 + double(-2.0) * rate_of_progress_1 + double(-2.0) * rate_of_progress_3 + double(-2.0) * rate_of_progress_4 + double(-2.0) * rate_of_progress_5 + double(-2.0) * rate_of_progress_6;
        net_production_rates[2] = double(1.0) * rate_of_progress_2;
        net_production_rates[3] = double(-1.0) * rate_of_progress_0 + double(-1.0) * rate_of_progress_1 + double(-1.0) * rate_of_progress_3 + double(-1.0) * rate_of_progress_4 + double(-1.0) * rate_of_progress_5 + double(-1.0) * rate_of_progress_6;
        //source_4 has no production term
        net_production_rates[5] = double(2.0) * rate_of_progress_0 + double(2.0) * rate_of_progress_1 + double(-1.0) * rate_of_progress_2 + double(2.0) * rate_of_progress_3 + double(2.0) * rate_of_progress_4 + double(2.0) * rate_of_progress_5 + double(2.0) * rate_of_progress_6;
        //source_6 has no production term
        //source_7 has no production term
        //source_8 has no production term
        //source_9 has no production term
        //source_10 has no production term
        //source_11 has no production term
        //source_12 has no production term
        //source_13 has no production term

        return net_production_rates;
    }
    
    SpeciesJacobian source_jacobian(const Species& species, const double& temperature)  
    {
        Species net_production_rates = {double(0)};
        SpeciesJacobian jacobian_net_production_rates = {double(0)};
        double inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature);
        double log_temperature = log_gen(temperature);
        double drate_of_progress_dspecies  = double(0);
        double equilibrium_constant  = double(0);
        Species drate_of_progress_dspecies_all_species  = {double(0)};
        
        Reactions gibbs_reactions = gibbs_reaction(log_temperature);
        Reactions dgibbs_reactions_dlog_temperature = dgibbs_reaction_dlog_temperature(log_temperature);

        
        double pressure_ = pressure(species, temperature);
        Species dpressure_dspecies_ = dpressure_dspecies(species, temperature); //unchecked
        
        double mixture_concentration = 
        multiply(pressure_,
                 inv_universal_gas_constant_temperature);
        
        Species dmixture_concentration_dspecies = Species{1}; // optimized (1/(RT))*(RT,...,RT)
        
            
        double forward_reaction_0 = call_forward_reaction_0(temperature, log_temperature);



        double forward_reaction_1 = call_forward_reaction_1(temperature, log_temperature);



        double forward_reaction_2 =  call_forward_reaction_2(species, temperature, log_temperature, mixture_concentration);

        Species dforward_reaction_2_dspecies = dcall_forward_reaction_2_dspecies(species,temperature,log_temperature,mixture_concentration);



        double forward_reaction_3 =  call_forward_reaction_3(temperature, pressure_);

        Species dforward_reaction_3_dspecies = scale_gen(dcall_forward_reaction_3_dpressure(temperature,pressure_), dpressure_dspecies_);


        double forward_reaction_4 =  call_forward_reaction_4(species, temperature, log_temperature, mixture_concentration);

        Species dforward_reaction_4_dspecies = dcall_forward_reaction_4_dspecies(species,temperature,log_temperature,mixture_concentration);



        double forward_reaction_5 =  call_forward_reaction_5(species, temperature, log_temperature, mixture_concentration);

        Species dforward_reaction_5_dspecies = dcall_forward_reaction_5_dspecies(species,temperature,log_temperature,mixture_concentration);



        double forward_reaction_6 =  call_forward_reaction_6(species, temperature, log_temperature, mixture_concentration);

        Species dforward_reaction_6_dspecies = dcall_forward_reaction_6_dspecies(species,temperature,log_temperature,mixture_concentration);



        equilibrium_constant = multiply(exp_gen(-gibbs_reactions[0]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
                //2 H2 + O2 <=> 2 H2O

        //drate_of_progress_temperature unused
                    //drate_of_progress_dspecies[0] = {double(0)};
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[1]) * species[3], forward_reaction_0); // [0][1] +
        jacobian_net_production_rates[1][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][1] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][1] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(pow_gen2(species[1]), forward_reaction_0); // [0][3] +
        jacobian_net_production_rates[1][3] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][3] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][3] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = -multiply(dpow_gen2_da(species[5]), divide(forward_reaction_0, equilibrium_constant));
        jacobian_net_production_rates[1][5] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][5] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][5] += 2.0*drate_of_progress_dspecies;
                //2 H2 + O2 => 2 H2O

                // rate_of_progress temperature derivative unused
                
        //drate_of_progress_dspecies[1] = {double(0)};
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[1]) * species[3], forward_reaction_1);
        jacobian_net_production_rates[1][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][1] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][1] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(pow_gen2(species[1]), forward_reaction_1);
        jacobian_net_production_rates[1][3] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][3] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][3] += 2.0*drate_of_progress_dspecies;
        equilibrium_constant = multiply(exp_gen(-gibbs_reactions[2]), pow_gen2(pressure_atmosphere() * inv_universal_gas_constant_temperature));
                //H2O + M <=> 2 H + O + M

        //drate_of_progress_temperature unused
            
        drate_of_progress_dspecies = forward_reaction_2;// 2 5 
        jacobian_net_production_rates[5][5] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[0][5] += 2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[2][5] += drate_of_progress_dspecies;
        drate_of_progress_dspecies = -multiply(dpow_gen2_da(species[0]) * species[2], divide(forward_reaction_2, equilibrium_constant));// 2 0
        jacobian_net_production_rates[5][0] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[0][0] += 2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[2][0] += drate_of_progress_dspecies;
        drate_of_progress_dspecies = -multiply(pow_gen2(species[0]), divide(forward_reaction_2, equilibrium_constant));// 2 2
        jacobian_net_production_rates[5][2] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[0][2] += 2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[2][2] += drate_of_progress_dspecies;


        drate_of_progress_dspecies_all_species = 
        scale_gen(species[5], dforward_reaction_2_dspecies) -
        scale_gen(divide(pow_gen2(species[0]) * species[2], 
                         equilibrium_constant), 
                  dforward_reaction_2_dspecies);
        jacobian_net_production_rates[5] = jacobian_net_production_rates[5] + scale_gen(-1.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[0] = jacobian_net_production_rates[0] + scale_gen(2.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[2] = jacobian_net_production_rates[2] + scale_gen(1.0, drate_of_progress_dspecies_all_species);

                                equilibrium_constant = multiply(exp_gen(-gibbs_reactions[3]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
                //2 H2 + O2 <=> 2 H2O

        //drate_of_progress_temperature unused
            
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[1]) * species[3], forward_reaction_3);//3 1
        jacobian_net_production_rates[1][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][1] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][1] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(pow_gen2(species[1]), forward_reaction_3);//3 3
        jacobian_net_production_rates[1][3] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][3] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][3] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = -multiply(dpow_gen2_da(species[5]), divide(forward_reaction_3, equilibrium_constant));// 3 5
        jacobian_net_production_rates[1][5] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][5] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][5] += 2.0*drate_of_progress_dspecies;


        drate_of_progress_dspecies_all_species = 
        scale_gen(pow_gen2(species[1]) * species[3], dforward_reaction_3_dspecies) -
        scale_gen(divide(pow_gen2(species[5]), 
                         equilibrium_constant), 
                  dforward_reaction_3_dspecies);
        jacobian_net_production_rates[1] = jacobian_net_production_rates[1] + scale_gen(-2.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[3] = jacobian_net_production_rates[3] + scale_gen(-1.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[5] = jacobian_net_production_rates[5] + scale_gen(2.0, drate_of_progress_dspecies_all_species);

                                equilibrium_constant = multiply(exp_gen(-gibbs_reactions[4]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
                //2 H2 + O2 (+M) <=> 2 H2O (+M)

        //drate_of_progress_temperature unused
            
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[1]) * species[3], forward_reaction_4);//4 1
        jacobian_net_production_rates[1][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][1] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][1] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(pow_gen2(species[1]), forward_reaction_4);//4 3
        jacobian_net_production_rates[1][3] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][3] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][3] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = -multiply(dpow_gen2_da(species[5]), divide(forward_reaction_4, equilibrium_constant));// 4 5
        jacobian_net_production_rates[1][5] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][5] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][5] += 2.0*drate_of_progress_dspecies;


        drate_of_progress_dspecies_all_species = 
        scale_gen(pow_gen2(species[1]) * species[3], dforward_reaction_4_dspecies) -
        scale_gen(divide(pow_gen2(species[5]), 
                         equilibrium_constant), 
                  dforward_reaction_4_dspecies);
        jacobian_net_production_rates[1] = jacobian_net_production_rates[1] + scale_gen(-2.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[3] = jacobian_net_production_rates[3] + scale_gen(-1.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[5] = jacobian_net_production_rates[5] + scale_gen(2.0, drate_of_progress_dspecies_all_species);

                                equilibrium_constant = multiply(exp_gen(-gibbs_reactions[5]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
                //2 H2 + O2 (+M) <=> 2 H2O (+M)

        //drate_of_progress_temperature unused
            
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[1]) * species[3], forward_reaction_5);//5 1
        jacobian_net_production_rates[1][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][1] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][1] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(pow_gen2(species[1]), forward_reaction_5);//5 3
        jacobian_net_production_rates[1][3] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][3] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][3] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = -multiply(dpow_gen2_da(species[5]), divide(forward_reaction_5, equilibrium_constant));// 5 5
        jacobian_net_production_rates[1][5] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][5] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][5] += 2.0*drate_of_progress_dspecies;


        drate_of_progress_dspecies_all_species = 
        scale_gen(pow_gen2(species[1]) * species[3], dforward_reaction_5_dspecies) -
        scale_gen(divide(pow_gen2(species[5]), 
                         equilibrium_constant), 
                  dforward_reaction_5_dspecies);
        jacobian_net_production_rates[1] = jacobian_net_production_rates[1] + scale_gen(-2.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[3] = jacobian_net_production_rates[3] + scale_gen(-1.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[5] = jacobian_net_production_rates[5] + scale_gen(2.0, drate_of_progress_dspecies_all_species);

                                equilibrium_constant = multiply(exp_gen(-gibbs_reactions[6]), inv_pressure_atmosphere() * universal_gas_constant() * temperature);
                //2 H2 + O2 (+M) <=> 2 H2O (+M)

        //drate_of_progress_temperature unused
            
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[1]) * species[3], forward_reaction_6);//6 1
        jacobian_net_production_rates[1][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][1] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][1] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(pow_gen2(species[1]), forward_reaction_6);//6 3
        jacobian_net_production_rates[1][3] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][3] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][3] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = -multiply(dpow_gen2_da(species[5]), divide(forward_reaction_6, equilibrium_constant));// 6 5
        jacobian_net_production_rates[1][5] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][5] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[5][5] += 2.0*drate_of_progress_dspecies;


        drate_of_progress_dspecies_all_species = 
        scale_gen(pow_gen2(species[1]) * species[3], dforward_reaction_6_dspecies) -
        scale_gen(divide(pow_gen2(species[5]), 
                         equilibrium_constant), 
                  dforward_reaction_6_dspecies);
        jacobian_net_production_rates[1] = jacobian_net_production_rates[1] + scale_gen(-2.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[3] = jacobian_net_production_rates[3] + scale_gen(-1.0, drate_of_progress_dspecies_all_species);
        jacobian_net_production_rates[5] = jacobian_net_production_rates[5] + scale_gen(2.0, drate_of_progress_dspecies_all_species);

                        

        /*
                //no temperature jacobian
        
            const int n_rates_of_progres_species_jacobian_0 = 1;
            static constexpr std::array<double, n_rates_of_progres_species_jacobian_0> coefficients_0 = {2.0};
            static constexpr std::array<int, n_rates_of_progres_species_jacobian_0> idx_0 = {2};
            
        for(int i = 0; i < n_species; i++)
        {
            double sum = 0;
            for(int j = 0; j < n_rates_of_progres_species_jacobian_0; j++)
            {
                sum += coefficients_0[j] * drate_of_progress_dspecies[idx_0[j]][i];
            }
            jacobian_net_production_rates[0][i] = sum;
        }
        */
            
        /*
                //no temperature jacobian
        
            const int n_rates_of_progres_species_jacobian_1 = 6;
            static constexpr std::array<double, n_rates_of_progres_species_jacobian_1> coefficients_1 = {-2.0,-2.0,-2.0,-2.0,-2.0,-2.0};
            static constexpr std::array<int, n_rates_of_progres_species_jacobian_1> idx_1 = {0,1,3,4,5,6};
            
        for(int i = 0; i < n_species; i++)
        {
            double sum = 0;
            for(int j = 0; j < n_rates_of_progres_species_jacobian_1; j++)
            {
                sum += coefficients_1[j] * drate_of_progress_dspecies[idx_1[j]][i];
            }
            jacobian_net_production_rates[1][i] = sum;
        }
        */
            
        /*
                //no temperature jacobian
        
            const int n_rates_of_progres_species_jacobian_2 = 1;
            static constexpr std::array<double, n_rates_of_progres_species_jacobian_2> coefficients_2 = {1.0};
            static constexpr std::array<int, n_rates_of_progres_species_jacobian_2> idx_2 = {2};
            
        for(int i = 0; i < n_species; i++)
        {
            double sum = 0;
            for(int j = 0; j < n_rates_of_progres_species_jacobian_2; j++)
            {
                sum += coefficients_2[j] * drate_of_progress_dspecies[idx_2[j]][i];
            }
            jacobian_net_production_rates[2][i] = sum;
        }
        */
            
        /*
                //no temperature jacobian
        
            const int n_rates_of_progres_species_jacobian_3 = 6;
            static constexpr std::array<double, n_rates_of_progres_species_jacobian_3> coefficients_3 = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
            static constexpr std::array<int, n_rates_of_progres_species_jacobian_3> idx_3 = {0,1,3,4,5,6};
            
        for(int i = 0; i < n_species; i++)
        {
            double sum = 0;
            for(int j = 0; j < n_rates_of_progres_species_jacobian_3; j++)
            {
                sum += coefficients_3[j] * drate_of_progress_dspecies[idx_3[j]][i];
            }
            jacobian_net_production_rates[3][i] = sum;
        }
        */
            
        //no species jacobian
            
        /*
                //no temperature jacobian
        
            const int n_rates_of_progres_species_jacobian_5 = 7;
            static constexpr std::array<double, n_rates_of_progres_species_jacobian_5> coefficients_5 = {2.0,2.0,-1.0,2.0,2.0,2.0,2.0};
            static constexpr std::array<int, n_rates_of_progres_species_jacobian_5> idx_5 = {0,1,2,3,4,5,6};
            
        for(int i = 0; i < n_species; i++)
        {
            double sum = 0;
            for(int j = 0; j < n_rates_of_progres_species_jacobian_5; j++)
            {
                sum += coefficients_5[j] * drate_of_progress_dspecies[idx_5[j]][i];
            }
            jacobian_net_production_rates[5][i] = sum;
        }
        */
            
        //no species jacobian
            
        //no species jacobian
            
        //no species jacobian
            
        //no species jacobian
            
        //no species jacobian
            
        //no species jacobian
            
        //no species jacobian
            
        //no species jacobian
            
        return jacobian_net_production_rates;
    }