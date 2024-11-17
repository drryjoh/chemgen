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
        auto source_threaded(const PointState& point_states) {const_option} 
        {{
            std::cout << "Beginning of source" <<std::endl;
            PointReactions point_reactions = std::make_unique<std::array<std::array<double, n_reactions>, n_points>>();
            PointReactions point_progress_rates = std::make_unique<std::array<std::array<double, n_reactions>, n_points>>();
            {index} chunk_size = 0;
            std::cout << "1" <<std::endl;
            std::unique_ptr<std::array<std::array<double, n_species>, n_points>> source;\n""".format(**vars(configuration)))
    
    def write_reaction_functions_pointer_list(self, file, reaction_calls, configuration):
        indentation = '        '
        reaction_call_list = ','.join([f'{indentation}reaction_{k}' for k in range(len(reaction_calls))])
        file.write("""
            {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}, {scalar_parameter}, {scalar_parameter}, {scalar_parameter}), n_reactions> reaction_functions = {{{reaction_call_list} }};
    """.format(**vars(configuration), reaction_call_list = reaction_call_list, n_reactions = len(reaction_calls)))

    def write_progress_rate_functions_pointer_list(self, file, reaction_calls, configuration):
        indentation = '        '
        progress_rates_call_list = ','.join([f'{indentation}progress_rate_{k}' for k in range(len(reaction_calls))])
        file.write("""
            {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}, {species_parameter}, {scalar_parameter}), n_reactions> progress_rate_functions = {{{progress_rates_call_list}}};
    """.format(**vars(configuration), progress_rates_call_list = progress_rates_call_list, n_reactions = len(reaction_calls)))

    def write_production_rate_functions_pointer_list(self, file, species_production_function_texts, configuration):
        indentation = '        '
        progress_rates_call_list = ','.join([f'{indentation}production_rate_{k}' for k in range(len(species_production_function_texts))])
        file.write("""
            {scalar_list}<{scalar} (*)({reactions_parameter}), n_species> production_rate_functions = {{{progress_rates_call_list}}};
    """.format(**vars(configuration), progress_rates_call_list = progress_rates_call_list, n_species = len(species_production_function_texts)))
    
    def write_forward_reaction_ttb_loop_with_timing(self, file, array, pointer_list, parameters, configuration):
        index = "{index}".format(**vars(configuration))
        scalar  = "{scalar}".format(**vars(configuration))
        species = "{species}".format(**vars(configuration))
        file.write(f"""
        // Measure serial execution time
        auto start_serial_{array} = std::chrono::high_resolution_clock::now();
        std::cout << "2" <<std::endl;
        for ({index} i = 0; i < n_reactions * n_points; ++i) 
        {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in

            // Static variables to retain values across iterations
            static Species species_;
            static double temperature_;
            static double log_temperature_;
            static double pressure_;
            static double mixture_concentration_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                log_temperature_ = log_gen(temperature_);
                pressure_ = pressure(species_, temperature_);
                mixture_concentration_ = mixture_concentration(species_);
            }}
            (*{array})[j][k] = {pointer_list}[k]({parameters});
        }}
        std::cout << "3" <<std::endl;

        auto end_serial_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time_{array} = end_serial_{array} - start_serial_{array};

        // Measure parallel execution time
        auto start_parallel_{array} = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(0, n_reactions * n_points, [&]({index} i) 
        {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            static Species species_;
            static double temperature_;
            static double log_temperature_;
            static double pressure_;
            static double mixture_concentration_;
            static Species gibbs_free_energy_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                log_temperature_ = log_gen(temperature_);
                pressure_ = pressure(species_, temperature_);
                mixture_concentration_ = mixture_concentration(species_);
            }}

            (*{array})[j][k] = {pointer_list}[k]({parameters});

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
            static Species species_;
            static double temperature_;
            static double log_temperature_;
            static double pressure_;
            static double mixture_concentration_;
            static Species gibbs_free_energy_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                log_temperature_ = log_gen(temperature_);
                pressure_ = pressure(species_, temperature_);
                mixture_concentration_ = mixture_concentration(species_);
            }}

            (*{array})[j][k] = {pointer_list}[k]({parameters});
            }}
        }});
        auto end_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array}_1 = end_parallel_{array}_1 - start_parallel_{array}_1;

        chunk_size = n_reactions * n_points / 128; ;  // Starting chunk size
        auto start_parallel_{array}_2 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_reactions * n_points, chunk_size), [&](const tbb::blocked_range<int>& r) {{
            for ({index} i = r.begin(); i < r.end(); ++i) 
            {{
            int j = i / n_reactions;  // Row index (point)
            int k = i % n_reactions;  // Column in (reaction)
            static Species species_;
            static double temperature_;
            static double log_temperature_;
            static double pressure_;
            static double mixture_concentration_;
            static Species gibbs_free_energy_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                log_temperature_ = log_gen(temperature_);
                pressure_ = pressure(species_, temperature_);
                mixture_concentration_ = mixture_concentration(species_);
            }}

            (*{array})[j][k] = {pointer_list}[k]({parameters});
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

    def write_progress_rates_ttb_loop_with_timing(self, file, array, pointer_list, parameters, configuration):
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

            // Static variables to retain values across iterations
            static Species species_;
            static double temperature_;
            static Species gibbs_free_energy_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                gibbs_free_energy_ = species_gibbs_energy_mole_specific(temperature_);
            }}
            auto point_reaction_ = (*point_reactions)[j][k];
            (*{array})[j][k] = {pointer_list}[k]({parameters});
        }}

        auto end_serial_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time_{array} = end_serial_{array} - start_serial_{array};

        // Measure parallel execution time
        auto start_parallel_{array} = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(0, n_reactions * n_points, [&]({index} i) 
        {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            static Species species_;
            static double temperature_;
            static Species gibbs_free_energy_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                gibbs_free_energy_ = species_gibbs_energy_mole_specific(temperature_);
            }}
            auto point_reaction_ = (*point_reactions)[j][k];
            (*{array})[j][k] = {pointer_list}[k]({parameters});

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
            static Species species_;
            static double temperature_;
            static Species gibbs_free_energy_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                gibbs_free_energy_ = species_gibbs_energy_mole_specific(temperature_);
            }}
            auto point_reaction_ = (*point_reactions)[j][k];
            (*{array})[j][k] = {pointer_list}[k]({parameters});
            }}
        }});
        auto end_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array}_1 = end_parallel_{array}_1 - start_parallel_{array}_1;

        chunk_size = n_reactions * n_points / 128; ;  // Starting chunk size
        auto start_parallel_{array}_2 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_reactions * n_points, chunk_size), [&](const tbb::blocked_range<int>& r) {{
            for ({index} i = r.begin(); i < r.end(); ++i) 
            {{
            int j = i / n_reactions;  // Row index (point)
            int k = i % n_reactions;  // Column in (reaction)
            static Species species_;
            static double temperature_;
            static Species gibbs_free_energy_;

            // Update only when k == 0
            if (k == 0) 
            {{
                species_ = species_from_chemical_state((*point_states)[j]);
                temperature_ = temperature_from_chemical_state((*point_states)[j]);
                gibbs_free_energy_ = species_gibbs_energy_mole_specific(temperature_);
            }}
            auto point_reaction_ = (*point_reactions)[j][k];
            (*{array})[j][k] = {pointer_list}[k]({parameters});
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

    def write_species_ttb_loop_with_timing(self, file, array, pointer_list, parameters, configuration):
        index = "{index}".format(**vars(configuration))
        file.write(f"""
        // Measure serial execution time
        //find me!
        auto start_serial_{array} = std::chrono::high_resolution_clock::now();
        for ({index} i = 0; i < n_species * n_points; ++i) {{
            int j = i / n_species;  // point
            int k = i % n_species;  // species
            auto progress_rates_ = (*point_progress_rates)[j];
            (*{array})[j][k] = {pointer_list}[k]({parameters});
        }}    
        auto end_serial_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time_{array} = end_serial_{array} - start_serial_{array};

        // Measure parallel execution time
        auto start_parallel_{array} = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(0, n_species * n_points, [&]({index} i) {{
            int j = i / n_species;  // point
            int k = i % n_species;  // species
            auto progress_rates_ = (*point_progress_rates)[j];
            (*{array})[j][k] = {pointer_list}[k]({parameters});
        }});
        auto end_parallel_{array} = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array} = end_parallel_{array} - start_parallel_{array};
        
        chunk_size = 20;  // Starting chunk size
        auto start_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_species * n_points, chunk_size), [&](const tbb::blocked_range<int>& r) {{
            for ({index} i = r.begin(); i < r.end(); ++i) {{
                int j = i / n_species;  // point
                int k = i % n_species;  // species
                auto progress_rates_ = (*point_progress_rates)[j];
                (*{array})[j][k] = {pointer_list}[k]({parameters});
            }}
        }});

        auto end_parallel_{array}_1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time_{array}_1 = end_parallel_{array}_1 - start_parallel_{array}_1;

        chunk_size = n_species * n_points/128;  // Starting chunk size
        auto start_parallel_{array}_2 = std::chrono::high_resolution_clock::now();
        tbb::parallel_for(tbb::blocked_range<int>(0, n_species * n_points, chunk_size), [&](const tbb::blocked_range<int>& r) {{
            for ({index} i = r.begin(); i < r.end(); ++i) {{
                int j = i / n_species;  // point
                int k = i % n_species;  // species
                auto progress_rates_ = (*point_progress_rates)[j];
                (*{array})[j][k] = {pointer_list}[k]({parameters});
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
        file.write(" return source;\n        }")

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

        self.write_forward_reaction_ttb_loop_with_timing(file, "point_reactions", "reaction_functions", "species_, temperature_, log_temperature_, pressure_, mixture_concentration_", configuration)
        self.write_progress_rates_ttb_loop_with_timing(file, "point_progress_rates", "progress_rate_functions", "species_, temperature_, gibbs_free_energy_, point_reaction_", configuration)
        self.write_species_ttb_loop_with_timing(file, "source", "production_rate_functions", "progress_rates_",configuration)

        self.write_end_of_function(file)
        headers.append('source.h')
    
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

    def write_reaction_ttb_loop(self, file, array, pointer_list, parameters, configuration):
        index = "{index}".format(**vars(configuration))
        file.write(f"""
        for ({index} i = 0; i < n_reactions * n_points; ++i) 
        {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            (*{array})[j][k] = {pointer_list}[k]({parameters});
        }}

        tbb::parallel_for(0, n_reactions * n_points, [&]({index} i) {{
            int j = i / n_reactions;  // Row index
            int k = i % n_reactions;  // Column in
            (*{array})[j][k] = {pointer_list}[k]({parameters});
        }});

        """)