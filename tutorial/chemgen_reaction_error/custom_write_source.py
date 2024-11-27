class SourceWriter:
    def write_start_of_source_function(self, file, configuration):
        file.write("""
        {device_option}
        {reactions_function} progress_rates({species_parameter} species, {scalar_parameter} temperature) {const_option} 
        {{
            {reactions} progress_rates = {{{scalar_cast}(0)}};
            {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
            {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);
            {scalar} log_temperature = log_gen(temperature);
            {scalar} pressure_ = pressure(species, temperature);
            {scalar} mixture_concentration = pressure_ * inv_universal_gas_constant_temperature;\n""".format(**vars(configuration)))
    
        file.write("\n")

    def write_reaction_calculations(self, file, reaction_calls, configuration):
        for reaction_index, reaction_call in enumerate(reaction_calls):
            file.write("            {scalar} forward_reaction_{reaction_index} = {reaction_call}".format(**vars(configuration), reaction_call=reaction_call, reaction_index = reaction_index))
    
    def write_progress_rates(self, file, progress_rates, is_reversible, equilibrium_constants, configuration):
        for i, progress_rate in enumerate(progress_rates):
            if is_reversible[i]:
                file.write("        {scalar} equilibrium_constant_{i} = {equilibrium_constant};\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration)))
            file.write(f"        {progress_rate}\n") 
        file.write("\n")

    def write_end_of_function(self, file, progress_rates):
        for k, progress_rate in enumerate(progress_rates):
            file.write(f"        progress_rates[{k}] = rate_of_progress_{k}; ")
        file.write("\n        return progress_rates;\n    }")

    def write_source(self, file, equilibrium_constants, 
                     reaction_calls,  progress_rates, is_reversible, species_production_on_fly_function_texts,
                     species_production_texts, headers, configuration): 
        self.write_start_of_source_function(file, configuration)
        self.write_reaction_calculations(file, reaction_calls, configuration)
        self.write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration)
        self.write_end_of_function(file, progress_rates)
        headers.append('source.h')
