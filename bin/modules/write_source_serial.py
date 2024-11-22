class SourceWriter:
    def write_start_of_source_function(self, file, configuration):
        file.write("""
        {device_option}
        {species_function} source({species_parameter} species, {scalar_parameter} temperature) {const_option} 
        {{
            {species} net_production_rates = {{{scalar_cast}(0)}};
            {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
            {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);
            {scalar} log_temperature = log_gen(temperature);
            {scalar} pressure_ = pressure(species, temperature);
            {scalar} mixture_concentration = pressure_ * inv_universal_gas_constant_temperature;\n""".format(**vars(configuration)))
    
    def write_progress_rates(self, file, progress_rates, is_reversible, equilibrium_constants, configuration):
        for i, progress_rate in enumerate(progress_rates):
            if is_reversible[i]:
                file.write("        {scalar} equilibrium_constant_{i} = {equilibrium_constant};\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration)))
            file.write(f"        {progress_rate}\n") 
        file.write("\n")
        
    def write_species_production(self, file, species_production_rates, configuration):
        for species_index, species_production in enumerate(species_production_rates):
            if species_production != '':
                file.write(f"        {configuration.source_element.format(i = species_index)} = {species_production};\n") 
            else:
                file.write(f"        //source_{species_index} has no production term\n")
        file.write("\n")

    def write_reaction_calculations(self, file, reaction_calls, configuration):
        for reaction_index, reaction_call in enumerate(reaction_calls):
            file.write("        {scalar} forward_reaction_{reaction_index} = {reaction_call}".format(**vars(configuration), reaction_call=reaction_call, reaction_index = reaction_index))


    def write_end_of_function(self, file):
        file.write("        return net_production_rates;\n    }")

    def write_source(self, file, equilibrium_constants, 
                     reaction_calls,  progress_rates, is_reversible, species_production_on_fly_function_texts,
                     species_production_texts, headers, configuration): 
        #
        self.write_start_of_source_function(file, configuration)
        self.write_reaction_calculations(file, reaction_calls, configuration)
        self.write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration)
        self.write_species_production(file, species_production_texts, configuration)
        headers.append('source.h')
        self.write_end_of_function(file)
