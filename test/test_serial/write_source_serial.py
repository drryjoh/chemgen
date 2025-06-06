class SourceWriter:
    def write_reaction_functions(self, file, reaction_calls, configuration):
        for reaction_index, reaction_call in enumerate(reaction_calls):
            file.write("{device_option}\n{scalar_function} reaction_{reaction_index}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} pressure_, {scalar_parameter} mixture_concentration) {const_option} {{ return {call} }}".format(**vars(configuration), reaction_index = reaction_index, call = reaction_call ))

    def write_progress_rate_functions(self, file, progress_rates, is_reversible, equilibrium_constants, configuration):
        for i, progress_rate in enumerate(progress_rates):
            if is_reversible[i]:
                file.write("        {device_option}\n{scalar_function} progress_rate_{i}({species_parameter} species, {scalar_parameter} temperature, {species_parameter} gibbs_free_energies, {scalar_parameter} forward_reaction_{i}) {const_option} {{{scalar} inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature); {scalar} equilibrium_constant_{i} =  {equilibrium_constant}; return {progress_rate}}}".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration), progress_rate = progress_rate.split("=")[1]))
            else:
                file.write("        {device_option}\n{scalar_function} progress_rate_{i}({species_parameter} species,{scalar_parameter} temperature, {species_parameter} gibbs_free_energies, {scalar_parameter} forward_reaction_{i}) {const_option} {{{scalar} inv_universal_gas_constant_temperature  = inv_gen(universal_gas_constant() * temperature); return {progress_rate}}}\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration), progress_rate = progress_rate.split("=")[1]))
        file.write("\n")

    def write_species_production_functions(self, file, species_production_rates, configuration):
        for species_index, species_production in enumerate(species_production_rates):
            if species_production != '':
                file.write("{device_option} {scalar_function} production_rate_{species_index}({reactions_parameter} progress_rates){{ return {species_production};}};\n".format(**vars(configuration), species_production = species_production, species_index = species_index)) 
            else:
                file.write(f"//source_{species_index} has no production term\n")
                file.write("{device_option} {scalar_function} production_rate_{species_index}({reactions_parameter} progress_rates){{ return {scalar_cast}(0);}};\n".format(**vars(configuration), species_production = species_production, species_index = species_index)) 
        file.write("\n")

    def write_start_of_source_function(self, file, configuration, fit_gibbs_reaction = True):
        file.write("""
        {device_option} {scalar_function} temperature_from_chemical_state( const ChemicalState& state) {const_option} {{ return state[0];}}
        
        {device_option} 
        {species_function} species_from_chemical_state( const ChemicalState& state) 
        {const_option} 
        {{ 
            Species species;
            std::copy(state.begin() + 1, state.end(), species.begin());
            return species;
        }}

        {device_option}
        PointSpecies source_serial(const PointState& point_states) {const_option} 
        {{
            PointSpecies  point_source = {{}};
            \n""".format(**vars(configuration)))
    
    def write_reaction_functions_pointer_list(self, file, reaction_calls, configuration):
        indentation = '        '
        reaction_call_list = ','.join([f'{indentation}reaction_{k}' for k in range(len(reaction_calls))])
        file.write("""
            {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}, {scalar_parameter}, {scalar_parameter}, {scalar_parameter}), n_reactions> reaction_functions = {{{reaction_call_list} }};
    """.format(**vars(configuration), reaction_call_list = reaction_call_list, n_reactions = len(reaction_calls)))

    def write_progress_rate_functions_pointer_list(self, file, reaction_calls, configuration):
        indentation = ''
        progress_rates_call_list = ','.join([f'{indentation}progress_rate_{k}' for k in range(len(reaction_calls))])
        file.write("""
            {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}, {species_parameter}, {scalar_parameter}), n_reactions> progress_rate_functions = {{{progress_rates_call_list}}};
    """.format(**vars(configuration), progress_rates_call_list = progress_rates_call_list, n_reactions = len(reaction_calls)))

    def write_production_rate_functions_pointer_list(self, file, species_production_function_texts, configuration):
        indentation = ''
        progress_rates_call_list = ','.join([f'{indentation}production_rate_{k}' for k in range(len(species_production_function_texts))])
        file.write("""
            {scalar_list}<{scalar} (*)({reactions_parameter}), n_species> production_rate_functions = {{{progress_rates_call_list}}};
    """.format(**vars(configuration), progress_rates_call_list = progress_rates_call_list, n_species = len(species_production_function_texts)))
    
    def write_point_loop(self, file, configuration):
        file.write("""
        
            for({index} i = 0; i < n_points; i++)
            {{
                auto state_ = point_states[i];
                auto temperature_ = temperature_from_chemical_state(state_);
                auto species_ = species_from_chemical_state(state_);
                auto gibbs_free_energies_ = species_gibbs_energy_mole_specific(temperature_);
                auto inv_universal_gas_constant_temperature_  = inv_gen(universal_gas_constant() * temperature_);
                auto log_temperature_ = log_gen(temperature_);
                auto pressure_ = pressure(species_, temperature_);
                auto mixture_concentration_ = pressure_ * inv_universal_gas_constant_temperature_;
                Reactions point_progress = {{}};
                for({index} j = 0; j < n_reactions; j++)
                {{
                    auto point_forward_reaction = reaction_functions[j](species_, temperature_, log_temperature_, pressure_, mixture_concentration_);
                    point_progress[j] = progress_rate_functions[j](species_, temperature_, gibbs_free_energies_, point_forward_reaction);
                }}
                for({index} j = 0; j < n_species; j++)
                {{
                    point_source[i][j] = production_rate_functions[j](point_progress);
                }}
            }} 
        """.format(**vars(configuration)))

    def write_end_of_function(self, file):
        file.write("        return point_source;\n    }")

    def write_source(self, file, equilibrium_constants, reaction_calls, 
                     progress_rates, is_reversible, species_production_on_fly_function_texts,
                     species_production_function_texts, headers, configuration, fit_gibbs_reaction = True): 
        self.write_reaction_functions(file, reaction_calls, configuration)
        self.write_progress_rate_functions(file, progress_rates, is_reversible, equilibrium_constants, configuration)
        self.write_species_production_functions(file, species_production_function_texts, configuration)

        self.write_start_of_source_function_threaded(file, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
        self.write_reaction_functions_pointer_list(file, reaction_calls, configuration)
        self.write_progress_rate_functions_pointer_list(file, reaction_calls, configuration)
        self.write_production_rate_functions_pointer_list(file, species_production_function_texts, configuration)
        self.write_point_loop(file, configuration)

        self.write_end_of_function(file)
        headers.append('source.h')