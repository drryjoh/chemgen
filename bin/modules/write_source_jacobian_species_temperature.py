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
        {scalar} drate_of_progress_dspecies  = {scalar_cast}(0);
        {scalar} equilibrium_constant  = {scalar_cast}(0);
        {scalar} dequilibrium_constant_dtemperature = {scalar_cast}(0);
        {species} drate_of_progress_dspecies_all_species  = {{{scalar_cast}(0)}};

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

        {species} dtemperature_dspecies_ = dtemperature_dspecies(species, temperature);
        
            \n""".format(**vars(configuration), gibbs = gibbs))
    def write_progress_jacobian_header(self, file, reaction_index, is_reversible, reactions_depend_on, configuration):
        pressure_dependency = ""
        if "pressure" in reactions_depend_on[reaction_index]:
            pressure_dependency = """\n                                               {scalar_parameter} dpressure_dtemperature_,
                                               {species_parameter} dpressure_dspecies_,\n""".format(**vars(configuration))
        else:
            pressure_dependency = ""

        if is_reversible[reaction_index]:
            file.write("""

void update_jacobian_reaction_{reaction_index}({jacobian}& jacobian_net_production_rates,
                                               {species_parameter} species,
                                               {scalar_parameter} temperature,
                                               {scalar_parameter} log_temperature,
                                               {scalar_parameter} mixture_concentration,
                                               {scalar_parameter} pressure_,{pressure_dependency}
                                               {species_parameter} dtemperature_dspecies_,
                                               {scalar_parameter} equilibrium_constant_{reaction_index},
                                               {scalar_parameter} dequilibrium_constant_{reaction_index}_dtemperature,
                                               {scalar_parameter} dlog_temperature_dtemperature)
{{
        """.format(reaction_index = reaction_index, **vars(configuration), pressure_dependency = pressure_dependency))
        else:
            file.write("""

void update_jacobian_reaction_{reaction_index}({jacobian}& jacobian_net_production_rates,
                                               {species_parameter} species,
                                               {scalar_parameter} temperature,
                                               {scalar_parameter} log_temperature,
                                               {scalar_parameter} mixture_concentration,
                                               {scalar_parameter} pressure_,{pressure_dependency}
                                               {species_parameter} dtemperature_dspecies_,
                                               {scalar_parameter} dlog_temperature_dtemperature)
{{
        """.format(reaction_index = reaction_index, **vars(configuration), pressure_dependency = pressure_dependency))
    
    def write_eq_and_derivatives(self, file, progress_rates, is_reversible, equilibrium_constants, dequilibrium_constants_dtemperature, configuration ):
        for i, progress_rate in enumerate(progress_rates):
            if is_reversible[i]:
                file.write("\n")
                file.write("        {scalar} equilibrium_constant_{i} = {equilibrium_constant};\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration)))
                file.write("        {scalar} dequilibrium_constant_{i}_dtemperature = {dequilibrium_constant};\n".format(i=i, dequilibrium_constant = dequilibrium_constants_dtemperature[i], **vars(configuration)))
                file.write("\n")

    def write_reaction_calculations_jacobian_i(self, file, reaction_calls, reactions_depend_on, reaction_index, configuration):
        reaction_call = reaction_calls[reaction_index]
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

    def write_progress_rates_jacobian(self, file, progress_rates, progress_rates_derivatives, reaction_calls, reactions_depend_on, is_reversible, configuration):

        for reaction_index, progress_rate in enumerate(progress_rates):
            self.write_progress_jacobian_header(file, reaction_index, is_reversible, reactions_depend_on, configuration)
            self.write_reaction_calculations_jacobian_i(file, reaction_calls, reactions_depend_on, reaction_index, configuration)
            file.write("{scalar} drate_of_progress_dspecies  = {scalar_cast}(0);\n".format(**vars(configuration)))
            file.write(f"        {progress_rate}\n")
            file.write(f"        {progress_rates_derivatives[reaction_index]}")
            file.write("}\n")

    def write_progress_rates_jacobian_calls(self, file, progress_rates, is_reversible, reactions_depend_on, configuration):
        for i, progress_rate in enumerate(progress_rates):
            pressure_dependency = ""
            if "pressure" in reactions_depend_on[i]:
                pressure_dependency = "dpressure_dtemperature_, dpressure_dspecies_,"
            if is_reversible[i]:
                file.write("""        update_jacobian_reaction_{reaction_index}(jacobian_net_production_rates, species, temperature, log_temperature, mixture_concentration, pressure_, {pressure_dependency}dtemperature_dspecies_, equilibrium_constant_{reaction_index}, dequilibrium_constant_{reaction_index}_dtemperature,dlog_temperature_dtemperature); \n""".format(reaction_index = i, pressure_dependency = pressure_dependency))
            else:
                file.write("""        update_jacobian_reaction_{reaction_index}(jacobian_net_production_rates, species, temperature, log_temperature, mixture_concentration, pressure_, {pressure_dependency}dtemperature_dspecies_, dlog_temperature_dtemperature); \n""".format(reaction_index = i, pressure_dependency = pressure_dependency))

        
    def write_species_production_jacobian(self, file, species_production_rates, configuration):
        print(species_production_rates)
        for species_index, species_production in enumerate(species_production_rates):
            if species_production != '':
                file.write(f"{species_production}") 
            else:
                file.write(f"//source_{species_index} has no production term\n")
        file.write("\n")

    def write_end_of_function_jacobian(self, file):
        file.write("        return jacobian_net_production_rates;\n    }")

    def write_source_species_temperature_derivative_header(self, file, reaction_index, is_reversible, reactions_depend_on, configuration):
        pressure_dependency = ""
        if "pressure" in reactions_depend_on[reaction_index]:
            pressure_dependency = """\n                                               {scalar_parameter} dpressure_dtemperature_,\n""".format(**vars(configuration))
        else:
            pressure_dependency = ""

        if is_reversible[reaction_index]:
            file.write("""
        void update_dsource_species_dtemperature_reaction_{reaction_index}({species}& dsource_species_dtemperature_,
                                                           {species_parameter} species,
                                                           {scalar_parameter} temperature,
                                                           {scalar_parameter} log_temperature,
                                                           {scalar_parameter} mixture_concentration,
                                                           {scalar_parameter} pressure_,{pressure_dependency}
                                                           {scalar_parameter} equilibrium_constant_{reaction_index},
                                                           {scalar_parameter} dequilibrium_constant_{reaction_index}_dtemperature,
                                                           {scalar_parameter} dlog_temperature_dtemperature)
{{
        """.format(reaction_index=reaction_index, **vars(configuration), pressure_dependency=pressure_dependency))
        else:
            file.write("""
        void update_dsource_species_dtemperature_reaction_{reaction_index}({species}& dsource_species_dtemperature_,
                                                           {species_parameter} species,
                                                           {scalar_parameter} temperature,
                                                           {scalar_parameter} log_temperature,
                                                           {scalar_parameter} mixture_concentration,
                                                           {scalar_parameter} pressure_,{pressure_dependency}
                                                           {scalar_parameter} dlog_temperature_dtemperature)
{{
        """.format(reaction_index=reaction_index, **vars(configuration), pressure_dependency=pressure_dependency))

    def write_source_species_temperature_derivative_i(self, file, reaction_calls, reactions_depend_on, reaction_index, configuration):
        reaction_call = reaction_calls[reaction_index]
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
            if dependent_variable == "log_temperature":
                file.write(f"        dforward_reaction_{reaction_index}_dtemperature += d{front}_dlog_temperature{back} * dlog_temperature_dtemperature;\n")
            if dependent_variable == "pressure":
                file.write(f"        dforward_reaction_{reaction_index}_dtemperature += d{front}_dpressure{back} * dpressure_dtemperature_;\n")
            file.write('\n')
        file.write('\n')

    def write_source_species_temperature_derivative(self, file, progress_rates, progress_rates_derivatives, reaction_calls, reactions_depend_on, is_reversible, configuration):
        file.write("""
        {device_option}
        {scalar_function}
        dspecies_internal_energy_mole_source_sum_dtemperature({species_parameter} species, {scalar_parameter} temperature, {species_parameter} dsource_species_dtemperature_) {const_option}
        {{
            return {sum}(molecular_weights() * multiply_chain(species_internal_energy_mass_specific(temperature),
                                                              dspecies_internal_energy_mass_specific_dtemperature(temperature),
                                                              source_species(species, temperature),
                                                              dsource_species_dtemperature_));
        }}

        {device_option}
        {scalar_function}
        dtemperature_source_dtemperature({scalar_parameter} temperature, {species_parameter} species, {scalar_parameter} dspecies_internal_energy_mole_source_sum_dtemperature_) {const_option}
        {{
            return
            -divide_chain(species_internal_energy_mole_source_sum(species, temperature),
                          dspecies_internal_energy_mole_source_sum_dtemperature_,
                          specific_heat_constant_volume_volume_specific(species, temperature),
                          dspecific_heat_constant_volume_volume_specific_dtemperature(species, temperature));
        }}

        {device_option}
        {species_function}
        dtemperature_source_dspecies({scalar_parameter} temperature, {species_parameter} species, {species_parameter} dspecies_internal_energy_mole_source_sum_dspecies_) {const_option}
        {{
            {scalar} alpha = species_internal_energy_mole_source_sum(species, temperature);
            {scalar} beta = specific_heat_constant_volume_volume_specific(species, temperature);

            return
            scale_gen(-inv_gen(beta), dspecies_internal_energy_mole_source_sum_dspecies_) + scale_gen(divide(alpha, pow2(beta)), dspecific_heat_constant_volume_volume_specific_dspecies(species, temperature));
        }}
        """.format(**vars(configuration)))

        for reaction_index, progress_rate in enumerate(progress_rates):
            self.write_source_species_temperature_derivative_header(file, reaction_index, is_reversible, reactions_depend_on, configuration)
            self.write_source_species_temperature_derivative_i(file, reaction_calls, reactions_depend_on, reaction_index, configuration)

            # Extract and remove temperature derivative stuff from progress_rate_derivatives
            progress_rate_temperature_derivative = []
            progress_rate_derivatives = progress_rates_derivatives[reaction_index]
            progress_rate_derivatives_split = progress_rate_derivatives.split(';')
            for i, line in reversed(list(enumerate(progress_rate_derivatives_split))): # reverse order so pop doesn't alter subsequent i; note that this creates a copy
                if "drate_of_progress_{reaction_index}_dtemperature".format(reaction_index=reaction_index) in line:
                    if "jacobian_net_production_rates" in line:
                        line = line.replace("jacobian_net_production_rates", "dsource_species_dtemperature_")
                        line = line.replace("[0]", "")
                        line = line.replace("] +=", "-1] +=")
                    progress_rate_temperature_derivative.append(line)
                    progress_rate_derivatives_split.pop(i)

            file.write(f"        {';'.join(progress_rate_temperature_derivative[::-1]) + ';'}")
            file.write("\n}\n")

            progress_rates_derivatives[reaction_index] = ';'.join(progress_rate_derivatives_split)

    def write_source_species_temperature_derivative_calls(self, file, progress_rates, is_reversible, reactions_depend_on, configuration):
        file.write("""\n        {species} dsource_species_dtemperature_ = {{{scalar_cast}(0)}};\n""".format(**vars(configuration)))
        for i, progress_rate in enumerate(progress_rates):
            pressure_dependency = ""
            if "pressure" in reactions_depend_on[i]:
                pressure_dependency = "dpressure_dtemperature_, dpressure_dspecies_,"
            if is_reversible[i]:
                file.write("""        update_dsource_species_dtemperature_reaction_{reaction_index}(dsource_species_dtemperature_, species, temperature, log_temperature, mixture_concentration, pressure_, {pressure_dependency}equilibrium_constant_{reaction_index}, dequilibrium_constant_{reaction_index}_dtemperature,dlog_temperature_dtemperature); \n""".format(reaction_index = i, pressure_dependency = pressure_dependency))
            else:
                file.write("""        update_dsource_species_dtemperature_reaction_{reaction_index}(dsource_species_dtemperature_, species, temperature, log_temperature, mixture_concentration, pressure_, {pressure_dependency}dlog_temperature_dtemperature); \n""".format(reaction_index = i, pressure_dependency = pressure_dependency))

        file.write("""
        jacobian_net_production_rates[0][0] = dtemperature_source_dtemperature(temperature, species, dspecies_internal_energy_mole_source_sum_dtemperature(species, temperature, dsource_species_dtemperature_));
""".format(**vars(configuration)))

        file.write("""
        {species} dspecies_internal_energy_mole_source_sum_dspecies_ = {{{scalar_cast}(0)}};
        {species} species_internal_energy_mole_ = molecular_weights() * species_internal_energy_mass_specific(temperature);

        for ({index} i = 0; i < n_species; i++)
        {{
            {species} jacobian_column;
            for ({index} j = 0; j < n_species; j++)
            {{
                jacobian_column[j] = jacobian_net_production_rates[j+1][i+1];
            }}
            dspecies_internal_energy_mole_source_sum_dspecies_[i] = dot(species_internal_energy_mole_, jacobian_column);
        }}

        {species} dtemperature_source_dspecies_ = dtemperature_source_dspecies(temperature, species, dspecies_internal_energy_mole_source_sum_dspecies_);
""".format(**vars(configuration)))

    def write_source_jacobian(self, file, equilibrium_constants, dequilibrium_constants_dtemperature, reactions_depend_on,
                     reaction_calls,  progress_rates, progress_rates_derivatives, is_reversible, species_production_on_fly_function_texts,
                     species_production_texts, species_production_jacobian_texts, headers, configuration, temperature_equation, fit_gibbs_reaction = True):
        if temperature_equation: self.write_source_species_temperature_derivative(file, progress_rates, progress_rates_derivatives, reaction_calls, reactions_depend_on, is_reversible, configuration)
        self.write_progress_rates_jacobian(file, progress_rates, progress_rates_derivatives, reaction_calls, reactions_depend_on, is_reversible, configuration)
        self.write_start_of_source_function_jacobian(file, configuration, fit_gibbs_reaction = fit_gibbs_reaction)
        self.write_eq_and_derivatives(file, progress_rates, is_reversible, equilibrium_constants, dequilibrium_constants_dtemperature, configuration)
        self.write_progress_rates_jacobian_calls(file, progress_rates, is_reversible, reactions_depend_on, configuration)
        if temperature_equation: self.write_source_species_temperature_derivative_calls(file, progress_rates, is_reversible, reactions_depend_on, configuration)
        self.write_end_of_function_jacobian(file)
        #headers.append('source.h')
