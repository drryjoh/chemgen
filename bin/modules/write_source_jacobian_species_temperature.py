class SourceJacobianWriter:
# Jacobian
    def write_start_of_source_function_jacobian(self, file, configuration, fit_gibbs_reaction = True):
        if fit_gibbs_reaction:
            gibbs = "{reactions} gibbs_reactions = gibbs_reaction(log_temperature);\n        {reactions} dgibbs_reactions_dlog_temperature = dgibbs_reaction_dlog_temperature(log_temperature);\n".format(**vars(configuration)) 
        else:
            gibbs = "{species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);".format(**vars(configuration)) 
        file.write("""
    {device_option}
    {jacobian_function} source_jacobian({species_parameter} species, {scalar_parameter} temperature) {const_option} 
    {{
        {species} net_production_rates = {{{scalar_cast}(0)}};
        {jacobian} jacobian_net_production_rates = {{{scalar_cast}(0)}};
        ReactionSpecies drate_of_progress_dspecies = {{{scalar_cast}(0)}};

        {scalar} inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature);
        {scalar} dinv_universal_gas_constant_temperature_dtemperature  = inv_chain(universal_gas_constant() * temperature, universal_gas_constant());
        
        {scalar} log_temperature = log_gen(temperature);
        {scalar} dlog_temperature_dtemperature = dlog_da(temperature);
        
        {gibbs}
        
        {scalar} pressure_ = pressure(species, temperature);
        {scalar} dpressure_dtemperature_ = dpressure_dtemperature(species, temperature); //unchecked
        {species} dpressure_dspecies_ = dpressure_dspecies(species, temperature); //unchecked
        
        {scalar} mixture_concentration = 
        multiply(pressure_,
                 inv_universal_gas_constant_temperature);
        {scalar} dmixture_concentration_dtemperature = 
        multiply_chain(pressure_, 
                       dpressure_dtemperature_,
                       inv_universal_gas_constant_temperature,
                       dinv_universal_gas_constant_temperature_dtemperature);
        
        {species} dmixture_concentration_dspecies = {species}{{1}}; // optimized (1/(RT))*(RT,...,RT)
        
            \n""".format(**vars(configuration), gibbs = gibbs))
    
    def write_progress_rates_jacobian(self, file, progress_rates, progress_rates_derivatives, is_reversible, equilibrium_constants, dequilibrium_constants_dtemperature, configuration):
        for i, progress_rate in enumerate(progress_rates):
            if is_reversible[i]:
                file.write("\n")
                file.write("        {scalar} equilibrium_constant_{i} = {equilibrium_constant};\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration)))
                file.write("        {scalar} dequilibrium_constant_{i}_dtemperature = {dequilibrium_constant};\n".format(i=i, dequilibrium_constant = dequilibrium_constants_dtemperature[i], **vars(configuration)))
                file.write("\n")
            file.write(f"        {progress_rate}\n")
            file.write(f"        {progress_rates_derivatives[i]}")
        file.write("\n")
        
    def write_species_production_jacobian(self, file, species_production_rates, configuration):
        for species_index, species_production in enumerate(species_production_rates):
            if species_production != '':
                file.write(f"{species_production}") 
            else:
                file.write(f"//source_{species_index} has no production term\n")
        file.write("\n")

    def write_reaction_calculations_jacobian(self, file, reaction_calls, reactions_depend_on, configuration):
        for reaction_index, reaction_call in enumerate(reaction_calls):
            file.write("        {scalar} forward_reaction_{reaction_index} = {reaction_call}".format(**vars(configuration), reaction_call=reaction_call, reaction_index = reaction_index))
            call  = f"{reaction_call}".replace(' ','')
            call_split = call.split('(')
            front = f"{call_split[0]}"
            back = f"({call_split[1]}".replace(';','').replace('\n','')

            #assure log_temperature and pressure come after species and temperature
            if "log_temperature" in reactions_depend_on[reaction_index]:
                reactions_depend_on[reaction_index].remove("log_temperature")
                reactions_depend_on[reaction_index].append("log_temperature")
            if "pressure" in reactions_depend_on[reaction_index]:
                reactions_depend_on[reaction_index].remove("pressure")
                reactions_depend_on[reaction_index].append("pressure")

            for dependent_variable in reactions_depend_on[reaction_index]:
                if dependent_variable == "temperature":
                    file.write(f"        {configuration.scalar} dforward_reaction_{reaction_index}_dtemperature = d{front}_dtemperature{back};\n")
                if dependent_variable == "species":
                    file.write(f"        {configuration.species} dforward_reaction_{reaction_index}_dspecies = d{front}_dspecies{back};\n")
                if dependent_variable == "log_temperature":
                    file.write(f"        dforward_reaction_{reaction_index}_dtemperature += d{front}_dlog_temperature{back} * dlog_temperature_dtemperature;\n")
                if dependent_variable == "pressure":
                    file.write(f"        dforward_reaction_{reaction_index}_dtemperature += d{front}_dpressure{back} * dpressure_dtemperature_;\n")
                    if "species" not in reactions_depend_on[reaction_index]:
                        file.write(f"        {configuration.species}   dforward_reaction_{reaction_index}_dspecies = scale_gen(d{front}_dpressure{back}, dpressure_dspecies_);\n")
                    else:
                        file.write(f"        dforward_reaction_{reaction_index}_dspecies += scale_gen(d{front}_dpressure{back}, dpressure_dspecies_);\n")
                file.write('\n')
            file.write('\n')

    def write_end_of_function_jacobian(self, file):
        file.write("        return jacobian_net_production_rates;\n    }")

    def write_source_jacobian(self, file, equilibrium_constants, dequilibrium_constants_dtemperature, reactions_depend_on,
                     reaction_calls,  progress_rates, progress_rates_derivatives, is_reversible, species_production_on_fly_function_texts,
                     species_production_texts, species_production_jacobian_texts, headers, configuration, fit_gibbs_reaction = True): 
        self.write_start_of_source_function_jacobian(file, configuration, fit_gibbs_reaction = fit_gibbs_reaction)
        self.write_reaction_calculations_jacobian(file, reaction_calls, reactions_depend_on, configuration)
        self.write_progress_rates_jacobian(file, progress_rates, progress_rates_derivatives, is_reversible, equilibrium_constants, dequilibrium_constants_dtemperature, configuration)
        self.write_species_production_jacobian(file, species_production_jacobian_texts, configuration)
        self.write_end_of_function_jacobian(file)
        #headers.append('source.h')
