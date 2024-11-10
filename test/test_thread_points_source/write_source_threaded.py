class SourceWriter:
    def write_reaction_functions(self, file, reaction_calls, configuration):
        for reaction_index, reaction_call in enumerate(reaction_calls):
            file.write("{device_option}\n{scalar_function} reaction_{reaction_index}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} pressure_, {scalar_parameter} mixture_concentration) {const_option} {{ return {call} }}".format(**vars(configuration), reaction_index = reaction_index, call = reaction_call ))

    def write_progress_rate_functions(self, file, progress_rates, is_reversible, equilibrium_constants, configuration):
        for i, progress_rate in enumerate(progress_rates):
            if is_reversible[i]:
                file.write("        {device_option}\n{scalar_function} progress_rate_{i}({species_parameter} species, {scalar_parameter} temperature, {species_parameter} gibbs_free_energies, {scalar_parameter} forward_reaction_{i}) {const_option} {{{scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature); {scalar} equilibrium_constant_{i} =  {equilibrium_constant}; return {progress_rate}}}".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration), progress_rate = progress_rate.split("=")[1]))
            else:
                file.write("        {device_option}\n{scalar_function} progress_rate_{i}({species_parameter} species,{scalar_parameter} temperature, {species_parameter} gibbs_free_energies, {scalar_parameter} forward_reaction_{i}) {const_option} {{{scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature); return {progress_rate}}}\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration), progress_rate = progress_rate.split("=")[1]))
        file.write("\n")

    def write_species_production_functions(self, file, species_production_rates, configuration):
        for species_index, species_production in enumerate(species_production_rates):
            if species_production != '':
                file.write("{device_option} {scalar_function} production_rate_{species_index}({reactions_parameter} progress_rates){{ return {species_production};}};\n".format(**vars(configuration), species_production = species_production, species_index = species_index)) 
            else:
                file.write(f"//source_{species_index} has no production term\n")
                file.write("{device_option} {scalar_function} production_rate_{species_index}({reactions_parameter} progress_rates){{ return {scalar_cast}(0);}};\n".format(**vars(configuration), species_production = species_production, species_index = species_index)) 
        file.write("\n")

    def write_start_of_source_function_threaded(self, file, configuration):
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
        {species_function} source_threaded(const PointState& point_states) {const_option} 
        {{
            PointState source = {{}};
            ChemicalState state = {{}};
            PointReactions point_reactions = {{}};
            PointReactions point_progress_rates = {{}};
            {species} net_production_rates = {{{scalar_cast}(0)}};
            {index} chunk_size = 0;
            std::array<std::array<double, n_species>, n_points> species, gibbs_free_energies;
            std::array<double, n_points> temperatures, pressures, log_temperatures, mixture_concentrations, inv_universal_gas_constant_temperatures;

            for({index} i = 0; i < n_points; i++)
            {{
                state = point_states[i];
                temperatures[i] = temperature_from_chemical_state(state);
                species[i] = species_from_chemical_state(state);
                gibbs_free_energies[i] = species_gibbs_energy_mole_specific(temperatures[i]);
                inv_universal_gas_constant_temperatures[i]  = inv(universal_gas_constant() * temperatures[i]);
                log_temperatures[i] = log_gen(temperatures[i]);
                pressures[i] = pressure(species[i], temperatures[i]);
                mixture_concentrations[i] = pressures[i] * inv_universal_gas_constant_temperatures[i];
            }} 

            {index} j = 0;
            {index} k = 0;\n""".format(**vars(configuration)))
    
    def write_reaction_functions_pointer_list(self, file, reaction_calls, configuration):
        indentation = '        '
        reaction_call_list = ','.join([f'\n{indentation}reaction_{k}' for k in range(len(reaction_calls))])
        file.write("""
            {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}, {scalar_parameter}, {scalar_parameter}, {scalar_parameter}), n_reactions> reaction_functions = {{{reaction_call_list} }};
    """.format(**vars(configuration), reaction_call_list = reaction_call_list, n_reactions = len(reaction_calls)))

    def write_progress_rate_functions_pointer_list(self, file, reaction_calls, configuration):
        indentation = '        '
        progress_rates_call_list = ','.join([f'\n{indentation}progress_rate_{k}' for k in range(len(reaction_calls))])
        file.write("""
            {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}, {species_parameter}, {scalar_parameter}), n_reactions> progress_rate_functions = {{{progress_rates_call_list}}};
    """.format(**vars(configuration), progress_rates_call_list = progress_rates_call_list, n_reactions = len(reaction_calls)))

    def write_production_rate_functions_pointer_list(self, file, species_production_function_texts, configuration):
        indentation = '        '
        progress_rates_call_list = ','.join([f'\n{indentation}production_rate_{k}' for k in range(len(species_production_function_texts))])
        file.write("""
            {scalar_list}<{scalar} (*)({reactions_parameter}), n_species> production_rate_functions = {{{progress_rates_call_list}}};
    """.format(**vars(configuration), progress_rates_call_list = progress_rates_call_list, n_species = len(species_production_function_texts)))

    def write_reaction_ttb_loop(self, file, array, pointer_list, parameters, configuration):
        index = "{index}".format(**vars(configuration))
        file.write(f"""
        for ({index} i = 0; i < n_reactions * n_points; ++i) 
        {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            {array}[j][k] = {pointer_list}[k]({parameters});
        }}

        tbb::parallel_for(0, n_reactions * n_points, [&]({index} i) {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            {array}[j][k] = {pointer_list}[k]({parameters});
        }});

        """)

    def write_reaction_ttb_loop_with_timing(self, file, array, pointer_list, parameters, configuration):
        index = "{index}".format(**vars(configuration))
        scalar  = "{scalar}".format(**vars(configuration))
        species = "{species}".format(**vars(configuration))
        file.write(f"""
        // Measure serial execution time
        auto start_serial_{array} = std::chrono::high_resolution_clock::now();
        for ({index} i = 0; i < n_reactions * n_points; ++i) 
        {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            auto species_ = species[j];
            auto temperature_ = temperatures[j];
            auto log_temperature_ = log_temperatures[j];
            auto pressure_ = pressures[j];
            auto mixture_concentration_ = mixture_concentrations[j];
            auto gibbs_free_energy_ = gibbs_free_energies[j];
            auto point_reaction_ = point_reactions[j][k];
            {array}[j][k] = {pointer_list}[k]({parameters});
        }}
        auto end_serial_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time_{array} = end_serial_{array} - start_serial_{array};

        // Measure parallel execution time
        auto start_parallel_{array} = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(0, n_reactions * n_points, [&]({index} i) 
        {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            auto species_ = species[j];
            auto temperature_ = temperatures[j];
            auto log_temperature_ = log_temperatures[j];
            auto pressure_ = pressures[j];
            auto mixture_concentration_ = mixture_concentrations[j];
            auto gibbs_free_energy_ = gibbs_free_energies[j];
            auto point_reaction_ = point_reactions[j][k];

            {array}[j][k] = {pointer_list}[k]({parameters});

        }});
        auto end_parallel_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array} = end_parallel_{array} - start_parallel_{array};

        chunk_size = 20;  // Starting chunk size
        auto start_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_reactions * n_points, chunk_size), [&](const tbb::blocked_range<int>& r) 
        {{
            for ({index} i = r.begin(); i < r.end(); ++i) 
            {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            auto species_ = species[j];
            auto log_temperature_ = log_temperatures[j];
            auto pressure_ = pressures[j];
            auto mixture_concentration_ = mixture_concentrations[j];
            auto temperature_ = temperatures[j];
            auto gibbs_free_energy_ = gibbs_free_energies[j];
            auto point_reaction_ = point_reactions[j][k];

            {array}[j][k] = {pointer_list}[k]({parameters});
            }}
        }});
        auto end_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array}_1 = end_parallel_{array}_1 - start_parallel_{array}_1;

        chunk_size = n_reactions * n_points / 16; ;  // Starting chunk size
        auto start_parallel_{array}_2 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_reactions * n_points, chunk_size), [&](const tbb::blocked_range<int>& r) {{
            for ({index} i = r.begin(); i < r.end(); ++i) 
            {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            auto species_ = species[j];
            auto temperature_ = temperatures[j];
            auto log_temperature_ = log_temperatures[j];
            auto pressure_ = pressures[j];
            auto mixture_concentration_ = mixture_concentrations[j];
            auto gibbs_free_energy_ = gibbs_free_energies[j];
            auto point_reaction_ = point_reactions[j][k];
            {array}[j][k] = {pointer_list}[k]({parameters});
            }}
        }});
        auto end_parallel_{array}_2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array}_2 = end_parallel_{array}_2 - start_parallel_{array}_2;

        // Output results
        std::cout << "Serial execution time for reactions loop {array}: " << serial_time_{array}.count() << " seconds"<<std::endl;
        std::cout << "Parallel execution time for reactions loop {array}: " << parallel_time_{array}.count() << " seconds"<<std::endl;
        std::cout << "Parallel execution time for reactions loop {array} chunk 20: " << parallel_time_{array}_1.count() << " seconds"<<std::endl;
        std::cout << "Parallel execution time for reactions loop {array} chunk 100: " << parallel_time_{array}_2.count() << " seconds"<<std::endl;

        """)

    def write_species_ttb_loop(self, file, array, pointer_list, parameters, configuration):
        index = "{index}".format(**vars(configuration))
        file.write(f"""
        
        for ({index} i = 0; i < n_species * n_points; ++i) {{
            {array}[i] = {pointer_list}[i]({parameters});
        }}

        tbb::parallel_for(0, n_species * n_points, [&]({index} i) {{
            {array}[i] = {pointer_list}[i]({parameters});
        }});

        """)

    def write_species_ttb_loop_with_timing(self, file, array, pointer_list, parameters, configuration):
        index = "{index}".format(**vars(configuration))
        file.write(f"""
        // Measure serial execution time
        auto start_serial_{array} = std::chrono::high_resolution_clock::now();
        for ({index} i = 0; i < n_species; ++i) {{
            {array}[i] = {pointer_list}[i]({parameters});
        }}    
        auto end_serial_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time_{array} = end_serial_{array} - start_serial_{array};

        // Measure parallel execution time
        auto start_parallel_{array} = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(0, n_species, [&]({index} i) {{
            {array}[i] = {pointer_list}[i]({parameters});
        }});
        auto end_parallel_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array} = end_parallel_{array} - start_parallel_{array};
        
        chunk_size = 20;  // Starting chunk size
        auto start_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_species, chunk_size), [&](const tbb::blocked_range<int>& r) {{
            for ({index} i = r.begin(); i < r.end(); ++i) {{
                {array}[i] = {pointer_list}[i]({parameters});
            }}
        }});
        auto end_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array}_1 = end_parallel_{array}_1 - start_parallel_{array}_1;

        chunk_size = 100;  // Starting chunk size
        auto start_parallel_{array}_2 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_species, chunk_size), [&](const tbb::blocked_range<int>& r) {{
            for ({index} i = r.begin(); i < r.end(); ++i) {{
                {array}[i] = {pointer_list}[i]({parameters});
            }}
        }});
        auto end_parallel_{array}_2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array}_2 = end_parallel_{array}_2 - start_parallel_{array}_2;
        // Output results
        std::cout << "Serial execution time for species_{array}: " << serial_time_{array}.count() << " seconds"<<std::endl;
        std::cout << "Parallel execution time for species_{array}: " << parallel_time_{array}.count() << " seconds"<<std::endl;
        std::cout << "Parallel execution time for reactions loop {array} chunk 20: " << parallel_time_{array}_1.count() << " seconds"<<std::endl;
        std::cout << "Parallel execution time for reactions loop {array} chunk 100: " << parallel_time_{array}_2.count() << " seconds"<<std::endl;

        """)

    def write_end_of_function(self, file):
        file.write("        return net_production_rates;\n    }")

    def write_source(self, file, equilibrium_constants, reaction_calls, 
                     progress_rates, is_reversible, 
                     species_production_function_texts, headers, configuration):
        self.write_reaction_functions(file, reaction_calls, configuration)
        self.write_progress_rate_functions(file, progress_rates, is_reversible, equilibrium_constants, configuration)
        self.write_species_production_functions(file, species_production_function_texts, configuration)

        self.write_start_of_source_function_threaded(file, configuration)

        self.write_reaction_functions_pointer_list(file, reaction_calls, configuration)
        self.write_progress_rate_functions_pointer_list(file, reaction_calls, configuration)
        self.write_production_rate_functions_pointer_list(file, species_production_function_texts, configuration)

        self.write_reaction_ttb_loop_with_timing(file, "point_reactions", "reaction_functions", "species_, temperature_, log_temperature_, pressure_, mixture_concentration_", configuration)
        self.write_reaction_ttb_loop_with_timing(file, "point_progress_rates", "progress_rate_functions", "species_, temperature_, gibbs_free_energy_, point_reaction_", configuration)
        #self.write_species_ttb_loop_with_timing(file, "net_production_rates", "production_rate_functions", "progress_rates",configuration)

        self.write_end_of_function(file)
        headers.append('source.h')