class SourceJacobianWriter:
# Jacobian
    def write_start_of_source_function_jacobian(self, file, configuration, fit_gibbs_reaction = True):
        if fit_gibbs_reaction:
            gibbs = "{reactions} gibbs_reactions = gibbs_reaction(log_temperature);\n        {reactions} dgibbs_reactions_dlog_temperature = dgibbs_reaction_dlog_temperature(log_temperature);\n".format(**vars(configuration)) 
        else:
            gibbs = "{species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);".format(**vars(configuration)) 
        file.write("""
    {device_option}
    {jacobian_function} source_jacobian_R({species_parameter} species, {scalar_parameter} temperature, {jacobian_parameter} J, {index} i) {const_option} 
    {{
        {species} net_production_rates = {{{scalar_cast}(0)}};
        {jacobian} jacobian_net_production_rates = J;
        {scalar} inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature);
        {scalar} log_temperature = log_gen(temperature);
        {scalar} drate_of_progress_dspecies  = {scalar_cast}(0);
        {scalar} equilibrium_constant  = {scalar_cast}(0);
        {species} drate_of_progress_dspecies_all_species  = {{{scalar_cast}(0)}};
        
        {gibbs}
        
        {scalar} pressure_ = pressure(species, temperature);
        {species} dpressure_dspecies_ = dpressure_dspecies(species, temperature); //unchecked
        
        {scalar} mixture_concentration = 
        multiply(pressure_,
                 inv_universal_gas_constant_temperature);
        
        {species} dmixture_concentration_dspecies = {species}{{1}}; // optimized (1/(RT))*(RT,...,RT)
        
            \n""".format(**vars(configuration), gibbs = gibbs))
    
    def write_progress_rates_jacobian(self, file, progress_rates, progress_rates_derivatives, is_reversible, equilibrium_constants, dequilibrium_constants_dtemperature, configuration):
        for i, progress_rate in enumerate(progress_rates):
            file.write(f"        if (i=={i}){{\n")
            if is_reversible[i]:
                file.write("        equilibrium_constant = {equilibrium_constant};\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration)))
            line = progress_rates_derivatives[i].replace("+=","-=")
            line = line.replace("] + s","] - s")
            file.write(f"        {line}\n        }}\n")
        file.write("\n")
        
    def write_species_production_jacobian(self, file, species_production_rates, configuration):
        #for species_index, species_production in enumerate(species_production_rates):
        #    if species_production != '':
        #        file.write(f"{species_production}") 
        #    else:
        #        file.write(f"//source_{species_index} has no production term\n")
        file.write("\n")

    def write_reaction_calculations_jacobian(self, file, reaction_calls, reactions_depend_on, configuration):
        for reaction_index, reaction_call in enumerate(reaction_calls):
            file.write("         {scalar} forward_reaction_{reaction_index} = {reaction_call}".format(**vars(configuration), reaction_call=reaction_call, reaction_index = reaction_index))
            call  = f"{reaction_call}".replace(' ','')
            call_split = call.split('(')
            front = f"{call_split[0]}"
            back = f"({call_split[1]}".replace(';','').replace('\n','')

            #assure log_temperature and pressure come after species and temperature
            if "log_temperature" in reactions_depend_on:
                reactions_depend_on.remove("log_temperature")
                reactions_depend_on.append("log_temperature")
            if "pressure" in reactions_depend_on:
                reactions_depend_on.remove("pressure")
                reactions_depend_on.append("pressure")

            dspecies_exist = False
            for dependent_variable in reactions_depend_on[reaction_index]:
                if dependent_variable == "species":
                    file.write(f"        {configuration.species} dforward_reaction_{reaction_index}_dspecies = d{front}_dspecies{back};\n")
                    dspecies_exist = True
                if dependent_variable == "pressure":
                    if dspecies_exist:
                        file.write(f"        dforward_reaction_{reaction_index}_dspecies += scale_gen(d{front}_dpressure{back}, dpressure_dspecies_);\n")
                    else:
                        file.write(f"        {configuration.species} dforward_reaction_{reaction_index}_dspecies = scale_gen(d{front}_dpressure{back}, dpressure_dspecies_);\n")
                file.write('\n')
            file.write('        \n')

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
