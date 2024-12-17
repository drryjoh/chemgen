class SourceWriter:
    def inline_functions(self, file, reaction_calls, progress_rates, is_reversible, equilibrium_constants, species_production_on_fly_function_texts, configuration):
        file.write("//parameters")
        file.write("//species, temperature, log_temperature, pressure, mixture_concentration, gibbs_free_energies, species_source")
        file.write("""    
    void compute_reaction ({species_parameter} species, 
                           {scalar_parameter} temperature, 
                           {scalar_parameter} log_temperature,
                           {scalar_parameter} pressure,
                           {scalar_parameter} mixture_concentration,
                           {species_parameter} gibbs_free_energies,
                           {species}& species_source)
    {{
""".format(**vars(configuration)))
        for reaction_index, reaction_call in enumerate(reaction_calls):
            if reaction_index == 0:
                scalar_declare  = "{scalar} ".format(**vars(configuration))
                file.write("       {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);\n".format(**vars(configuration)))
            else:
                scalar_declare  = ""

            progress_rate = progress_rates[reaction_index]
            species_production = species_production_on_fly_function_texts[reaction_index]
            file.write(f"//inline_{reaction_index}\n")
            
            file.write("        {scalar_declare}forward_reaction = {call}".format(**vars(configuration), reaction_index = reaction_index, call = reaction_call.replace("pressure_","pressure"), scalar_declare = scalar_declare ))
            if is_reversible[reaction_index]:
                file.write("        {scalar_declare}equilibrium_constant =  {equilibrium_constant};\n        {scalar_declare}progress_rate = {progress_rate}".format(reaction_index = reaction_index, equilibrium_constant = equilibrium_constants[reaction_index], **vars(configuration), progress_rate = progress_rate.split("=")[1].replace(f'_{reaction_index}',''), scalar_declare = scalar_declare))
            else:
                file.write("        {scalar_declare}progress_rate = {progress_rate}\n".format(reaction_index=reaction_index, equilibrium_constant = equilibrium_constants[reaction_index], **vars(configuration), progress_rate = progress_rate.split("=")[1].replace(f"_{reaction_index}",""), scalar_declare = scalar_declare))
            file.write("""
        {species_production}
""".format(**vars(configuration), species_production = species_production.replace('\n','\n        '), reaction_index = reaction_index)) 
        file.write("}\n")

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
        auto source_decomposed(int argc, char* argv[], const PointState& point_states) {const_option} 
        {{
            PointSpecies point_source = std::make_unique<std::array<std::array<double, n_species>, n_points>>();
            MPI_Init(&argc, &argv);

            {index} rank, n_ranks;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

            // Determine the number of points each rank will handle
            {index} points_per_rank = n_points / n_ranks;
            {index} remainder = n_points % n_ranks;

            // Calculate start and end index for each rank
            {index} start = rank * points_per_rank + std::min(rank, remainder);
            {index} end = start + points_per_rank + (rank < remainder ? 1 : 0);
            {index} local_n_points = points_per_rank + (rank < remainder ? 1 : 0);

            std::vector<{scalar_list}<{scalar}, n_species>> point_source_local(local_n_points, {scalar_list}<double, n_species>{{0.0}});


            \n""".format(**vars(configuration)))
    
    def write_point_loop(self, file, configuration):
        file.write("""
        auto start_serial = std::chrono::high_resolution_clock::now();
        {index} max_threads = tbb::this_task_arena::max_concurrency();
        if (rank == 0) 
        {{
            std::cout << "Number of threads available: " << max_threads << std::endl;
        }}
        {index} chunk_size = (end - start)/(max_threads);  // Starting chunk size
        tbb::parallel_for(tbb::blocked_range<int>(start, end, chunk_size), [&](const tbb::blocked_range<int>& r) 
        {{
            for ({index} i = r.begin(); i < r.end(); ++i) 
            {{
                auto state_ = (*point_states)[i];
                auto temperature_ = temperature_from_chemical_state(state_);
                auto species_ = species_from_chemical_state(state_);
                auto gibbs_free_energies_ = species_gibbs_energy_mole_specific(temperature_);
                auto inv_universal_gas_constant_temperature_  = inv(universal_gas_constant() * temperature_);
                auto log_temperature_ = log_gen(temperature_);
                auto pressure_ = pressure(species_, temperature_);
                auto mixture_concentration_ = pressure_ * inv_universal_gas_constant_temperature_;
                compute_reaction(species_, 
                                    temperature_, 
                                    log_temperature_, 
                                    pressure_,
                                    mixture_concentration_, 
                                    gibbs_free_energies_,
                                    point_source_local[i - start]);
            }}
        }}, tbb::auto_partitioner{{}} );
        
        auto end_serial = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time = end_serial - start_serial; 
        """.format(**vars(configuration)))

    def write_end_of_function(self, file, configuration):
        mpi_scalar_type = "MPI_{}".format("{scalar}".format(**vars(configuration)).upper())
        file.write("""

            if (rank == 0) {{
                // Copy local results directly for rank 0
                for (int i = 0; i < points_per_rank + (remainder > 0 ? 1 : 0); ++i) {{
                    (*point_source)[i] = point_source_local[i];
                }}
                // Receive data from other ranks
                for (int r = 1; r < n_ranks; ++r) {{
                    int recv_start = r * points_per_rank + std::min(r, remainder);
                    int recv_end = recv_start + points_per_rank + (r < remainder ? 1 : 0);
                    // Use reinterpret_cast to access nested arrays correctly
                    MPI_Recv(reinterpret_cast<double*>(&(*point_source)[recv_start][0]),
                            (recv_end - recv_start) * n_species, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }}
            }} else {{
                // Send local results to the root process
                MPI_Send(reinterpret_cast<const double*>(&point_source_local[0][0]),
                        local_n_points * n_species, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }}


            if (rank == 0) {{
                std::cout << "Computation completed and gathered on root process."<<std::endl;
                std::cout << "Rank 1 execution time: " << serial_time.count() << " seconds"<<std::endl;
            }}

            MPI_Finalize();
            return point_source;\n    }}""".format(**vars(configuration), mpi_scalar_type = mpi_scalar_type))

    def write_source(self, file, equilibrium_constants, reaction_calls, 
                     progress_rates, is_reversible, species_production_on_fly_function_texts,
                     species_production_function_texts, headers, configuration, fit_gibbs_reaction = True): 
        self.inline_functions(file, reaction_calls, progress_rates, 
                              is_reversible, equilibrium_constants, species_production_on_fly_function_texts, configuration)

        self.write_start_of_source_function_threaded(file, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
        self.write_point_loop(file, configuration)

        self.write_end_of_function(file, configuration)
        headers.append('source.h')