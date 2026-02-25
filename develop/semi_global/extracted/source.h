
        
        Species source_species(const Species& species, const double& temperature) 
        {
            Species net_production_rates = {double(0)};
            double inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature);
            double log_temperature = log_gen(temperature);
            Reactions gibbs_reactions = gibbs_reaction(log_temperature);
            double pressure_ = pressure(species, temperature);
            double mixture_concentration = pressure_ * inv_universal_gas_constant_temperature;
			
			//add the species^power here
			double forward_reaction_0 = call_forward_reaction_0(temperature, log_temperature);
			double forward_reaction_1 = call_forward_reaction_1(temperature, log_temperature);
			double forward_reaction_2 = call_forward_reaction_2(temperature, log_temperature);
			
			//I can't remeber: for semi global does the species^power above take care of the the stoichemetric powers below?
			double rate_of_progress_0 = multiply(pow_gen2(species[1]) * pow_gen7(species[0]), forward_reaction_0);
			double rate_of_progress_1 = multiply(pow_gen2(species[2]) * species[0], forward_reaction_1);
			double rate_of_progress_2 = multiply(pow_gen2(species[3]), forward_reaction_2);

			net_production_rates[0] = double(-7.0) * rate_of_progress_0 + double(-1.0) * rate_of_progress_1 + double(1.0) * rate_of_progress_2;
			net_production_rates[1] = double(-2.0) * rate_of_progress_0;
			net_production_rates[2] = double(8.0) * rate_of_progress_0 + double(-2.0) * rate_of_progress_1 + double(2.0) * rate_of_progress_2;
			net_production_rates[3] = double(2.0) * rate_of_progress_1 + double(-2.0) * rate_of_progress_2;
			net_production_rates[4] = double(6.0) * rate_of_progress_0;

			return net_production_rates;
		}
		
        double source_energy(const Species& species, const double& temperature) 
		{
			return double(0);
		}
        
		//source is (rho*u, C_1, ..., C_{n_s})
        ChemicalState source(const Species& species, const double& temperature) 
        {
            return set_chemical_state(source_energy(species, temperature), source_species(species, temperature));
        }

//no need for analytical jacobian now

#if 0            
        void update_dsource_species_dtemperature_reaction_0(Species& dsource_species_dtemperature_,
                                                           const Species& species,
                                                           const double& temperature,
                                                           const double& log_temperature,
                                                           const double& mixture_concentration,
                                                           const double& pressure_,
                                                           const double& dlog_temperature_dtemperature)
{
                double forward_reaction_0 = call_forward_reaction_0(temperature, log_temperature);
        double dforward_reaction_0_dtemperature = dcall_forward_reaction_0_dtemperature(temperature,log_temperature);

        dforward_reaction_0_dtemperature += dcall_forward_reaction_0_dlog_temperature(temperature,log_temperature) * dlog_temperature_dtemperature;


        

        // Reaction 0: 2 C4H6 + 7 O2 => 8 CO + 6 H2O

        double drate_of_progress_0_dtemperature =  multiply(pow_gen2(species[1]) * pow_gen7(species[0]), dforward_reaction_0_dtemperature);
        dsource_species_dtemperature_[1] += -2.0*drate_of_progress_0_dtemperature;
        dsource_species_dtemperature_[0] += -7.0*drate_of_progress_0_dtemperature;
        dsource_species_dtemperature_[2] += 8.0*drate_of_progress_0_dtemperature;
        dsource_species_dtemperature_[4] += 6.0*drate_of_progress_0_dtemperature;
}

        void update_dsource_species_dtemperature_reaction_1(Species& dsource_species_dtemperature_,
                                                           const Species& species,
                                                           const double& temperature,
                                                           const double& log_temperature,
                                                           const double& mixture_concentration,
                                                           const double& pressure_,
                                                           const double& dlog_temperature_dtemperature)
{
                double forward_reaction_1 = call_forward_reaction_1(temperature, log_temperature);
        double dforward_reaction_1_dtemperature = dcall_forward_reaction_1_dtemperature(temperature,log_temperature);

        dforward_reaction_1_dtemperature += dcall_forward_reaction_1_dlog_temperature(temperature,log_temperature) * dlog_temperature_dtemperature;


        

        // Reaction 1: 2 CO + O2 => 2 CO2

        double drate_of_progress_1_dtemperature =  multiply(pow_gen2(species[2]) * species[0], dforward_reaction_1_dtemperature);
        dsource_species_dtemperature_[2] += -2.0*drate_of_progress_1_dtemperature;
        dsource_species_dtemperature_[0] += -drate_of_progress_1_dtemperature;
        dsource_species_dtemperature_[3] += 2.0*drate_of_progress_1_dtemperature;
}

        void update_dsource_species_dtemperature_reaction_2(Species& dsource_species_dtemperature_,
                                                           const Species& species,
                                                           const double& temperature,
                                                           const double& log_temperature,
                                                           const double& mixture_concentration,
                                                           const double& pressure_,
                                                           const double& dlog_temperature_dtemperature)
{
                double forward_reaction_2 = call_forward_reaction_2(temperature, log_temperature);
        double dforward_reaction_2_dtemperature = dcall_forward_reaction_2_dtemperature(temperature,log_temperature);

        dforward_reaction_2_dtemperature += dcall_forward_reaction_2_dlog_temperature(temperature,log_temperature) * dlog_temperature_dtemperature;


        

        // Reaction 2: 2 CO2 => 2 CO + O2

        double drate_of_progress_2_dtemperature =  multiply(pow_gen2(species[3]), dforward_reaction_2_dtemperature);
        dsource_species_dtemperature_[3] += -2.0*drate_of_progress_2_dtemperature;
        dsource_species_dtemperature_[2] += 2.0*drate_of_progress_2_dtemperature;
        dsource_species_dtemperature_[0] += drate_of_progress_2_dtemperature;
}


void update_jacobian_reaction_0(SpeciesJacobian& jacobian_net_production_rates,
                                               const Species& species,
                                               const double& temperature,
                                               const double& log_temperature,
                                               const double& mixture_concentration,
                                               const double& pressure_,
                                               const Species& dtemperature_dspecies_,
                                               const double& dlog_temperature_dtemperature)
{
                double forward_reaction_0 = call_forward_reaction_0(temperature, log_temperature);
        double dforward_reaction_0_dtemperature = dcall_forward_reaction_0_dtemperature(temperature,log_temperature);

        dforward_reaction_0_dtemperature += dcall_forward_reaction_0_dlog_temperature(temperature,log_temperature) * dlog_temperature_dtemperature;


double drate_of_progress_dspecies  = double(0);
        double rate_of_progress_0 = multiply(pow_gen2(species[1]) * pow_gen7(species[0]), forward_reaction_0);
        


        //drate_of_progress_dspecies[0] = {double(0)};
        drate_of_progress_dspecies = multiply(dpow_gen7_da(species[0]) * pow_gen2(species[1]), forward_reaction_0);
        jacobian_net_production_rates[2][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[1][1] += -7.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][1] += 8.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[5][1] += 6.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[1]) * pow_gen7(species[0]), forward_reaction_0);
        jacobian_net_production_rates[2][2] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[1][2] += -7.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][2] += 8.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[5][2] += 6.0*drate_of_progress_dspecies;
}


void update_jacobian_reaction_1(SpeciesJacobian& jacobian_net_production_rates,
                                               const Species& species,
                                               const double& temperature,
                                               const double& log_temperature,
                                               const double& mixture_concentration,
                                               const double& pressure_,
                                               const Species& dtemperature_dspecies_,
                                               const double& dlog_temperature_dtemperature)
{
                double forward_reaction_1 = call_forward_reaction_1(temperature, log_temperature);
        double dforward_reaction_1_dtemperature = dcall_forward_reaction_1_dtemperature(temperature,log_temperature);

        dforward_reaction_1_dtemperature += dcall_forward_reaction_1_dlog_temperature(temperature,log_temperature) * dlog_temperature_dtemperature;


double drate_of_progress_dspecies  = double(0);
        double rate_of_progress_1 = multiply(pow_gen2(species[2]) * species[0], forward_reaction_1);
        


        //drate_of_progress_dspecies[1] = {double(0)};
        drate_of_progress_dspecies = multiply(pow_gen2(species[2]), forward_reaction_1);
        jacobian_net_production_rates[3][1] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[1][1] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[4][1] += 2.0*drate_of_progress_dspecies;
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[2]) * species[0], forward_reaction_1);
        jacobian_net_production_rates[3][3] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[1][3] += -drate_of_progress_dspecies;
        jacobian_net_production_rates[4][3] += 2.0*drate_of_progress_dspecies;
}


void update_jacobian_reaction_2(SpeciesJacobian& jacobian_net_production_rates,
                                               const Species& species,
                                               const double& temperature,
                                               const double& log_temperature,
                                               const double& mixture_concentration,
                                               const double& pressure_,
                                               const Species& dtemperature_dspecies_,
                                               const double& dlog_temperature_dtemperature)
{
                double forward_reaction_2 = call_forward_reaction_2(temperature, log_temperature);
        double dforward_reaction_2_dtemperature = dcall_forward_reaction_2_dtemperature(temperature,log_temperature);

        dforward_reaction_2_dtemperature += dcall_forward_reaction_2_dlog_temperature(temperature,log_temperature) * dlog_temperature_dtemperature;


double drate_of_progress_dspecies  = double(0);
        double rate_of_progress_2 = multiply(pow_gen2(species[3]), forward_reaction_2);
        


        //drate_of_progress_dspecies[2] = {double(0)};
        drate_of_progress_dspecies = multiply(dpow_gen2_da(species[3]), forward_reaction_2);
        jacobian_net_production_rates[4][4] += -2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[3][4] += 2.0*drate_of_progress_dspecies;
        jacobian_net_production_rates[1][4] += drate_of_progress_dspecies;
}

    
    SpeciesJacobian source_jacobian(const Species& species, const double& temperature)  
    {
        Species net_production_rates = {double(0)};
        SpeciesJacobian jacobian_net_production_rates = {double(0)};
        double drate_of_progress_dspecies  = double(0);
        double equilibrium_constant  = double(0);
        double dequilibrium_constant_dtemperature = double(0);
        Species drate_of_progress_dspecies_all_species  = {double(0)};

        double inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature);
        double dinv_universal_gas_constant_temperature_dtemperature  = inv_chain(universal_gas_constant() * temperature, universal_gas_constant());
        
        double log_temperature = log_gen(temperature);
        double dlog_temperature_dtemperature = dlog_da(temperature);
        
        Reactions gibbs_reactions = gibbs_reaction(log_temperature);
        Reactions dgibbs_reactions_dlog_temperature = dgibbs_reaction_dlog_temperature(log_temperature);

        
        double pressure_ = pressure(species, temperature);
        double dpressure_dtemperature_ = dpressure_dtemperature(species, temperature); //unchecked
        Species dpressure_dspecies_ = dpressure_dspecies(species, temperature); //unchecked
        
        double mixture_concentration = 
        multiply(pressure_,
                 inv_universal_gas_constant_temperature);
        double dmixture_concentration_dtemperature = 
        multiply_chain(pressure_, 
                       dpressure_dtemperature_,
                       inv_universal_gas_constant_temperature,
                       dinv_universal_gas_constant_temperature_dtemperature);
        
        Species dmixture_concentration_dspecies = Species{1}; // optimized (1/(RT))*(RT,...,RT)

        Species dtemperature_dspecies_ = dtemperature_dspecies(species, temperature);
        
            
        update_jacobian_reaction_0(jacobian_net_production_rates, species, temperature, log_temperature, mixture_concentration, pressure_, dtemperature_dspecies_, dlog_temperature_dtemperature); 
        update_jacobian_reaction_1(jacobian_net_production_rates, species, temperature, log_temperature, mixture_concentration, pressure_, dtemperature_dspecies_, dlog_temperature_dtemperature); 
        update_jacobian_reaction_2(jacobian_net_production_rates, species, temperature, log_temperature, mixture_concentration, pressure_, dtemperature_dspecies_, dlog_temperature_dtemperature); 

        Species dsource_species_dtemperature_ = {double(0)};
        update_dsource_species_dtemperature_reaction_0(dsource_species_dtemperature_, species, temperature, log_temperature, mixture_concentration, pressure_, dlog_temperature_dtemperature); 
        update_dsource_species_dtemperature_reaction_1(dsource_species_dtemperature_, species, temperature, log_temperature, mixture_concentration, pressure_, dlog_temperature_dtemperature); 
        update_dsource_species_dtemperature_reaction_2(dsource_species_dtemperature_, species, temperature, log_temperature, mixture_concentration, pressure_, dlog_temperature_dtemperature); 

        for (int i = 0; i < n_species; i++)
        {
            for (int j = 0; j < n_species; j++)
            {
                // temperature dependence
                jacobian_net_production_rates[i+1][j+1] += scale_gen(dsource_species_dtemperature_[i], dtemperature_dspecies_[j]);
            }
        }

        return jacobian_net_production_rates;
    }
#endif
