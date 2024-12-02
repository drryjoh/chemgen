class SourceWriter:
    def write_start_of_source_function(self, file, configuration):
        file.write("""
        {device_option}
        auto reactions_I_want({species_parameter} species, {scalar_parameter} temperature) {const_option} 
        {{
            ReactionSubset reaction_subset = {{{scalar_cast}(0)}};
            {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
            {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);
            {scalar} log_temperature = log_gen(temperature);
            {scalar} pressure_ = pressure(species, temperature);
            {scalar} mixture_concentration = pressure_ * inv_universal_gas_constant_temperature;\n""".format(**vars(configuration)))
    
        file.write("\n")

    def write_reaction_calculations(self, file, reaction_calls, configuration):
        reactions_of_interest = list(set([157, 183, 430 ,634, 945]))
        for k, reaction_index in enumerate(reactions_of_interest):
                file.write("            reaction_subset[{k}] = {reaction_call}".format(**vars(configuration), reaction_call=reaction_calls[reaction_index], k=k))
        file.write("return reaction_subset;\n}")

    def write_source(self, file, equilibrium_constants, 
                     reaction_calls,  progress_rates, is_reversible, species_production_on_fly_function_texts,
                     species_production_texts, headers, configuration): 
        self.write_start_of_source_function(file, configuration)
        self.write_reaction_calculations(file, reaction_calls, configuration)
        headers.append('source.h')
