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

    def write_species_production_functions_on_the_fly(self, file, species_production_on_fly_function_texts, configuration):
        for reaction_index, species_production in enumerate(species_production_on_fly_function_texts):
            file.write("""
{device_option} 
void production_rate_{reaction_index}({scalar_parameter} progress_rate, {species}& species_source)
{{
{species_production}
}};
""".format(**vars(configuration), species_production = species_production, reaction_index = reaction_index)) 

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

    def write_accrue_source_pointer_list(self, file, species_production_on_fly_function_texts, configuration):
        indentation = ''
        progress_rates_call_list = ','.join([f'{indentation}production_rate_{k}' for k in range(len(species_production_on_fly_function_texts))])
        file.write("""
            {scalar_list}<void (*)({scalar_parameter}, {species}&), n_reactions> accrue_source_on_fly = {{{progress_rates_call_list}}};
    """.format(**vars(configuration), progress_rates_call_list = progress_rates_call_list))
    
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
                auto inv_universal_gas_constant_temperature_  = inv_gen(universal_gas_constant() * temperature_);
                auto log_temperature_ = log_gen(temperature_);
                auto pressure_ = pressure(species_, temperature_);
                auto mixture_concentration_ = pressure_ * inv_universal_gas_constant_temperature_;
                for({index} j = 0; j < n_reactions; j++)
                {{
                    auto point_forward_reaction = reaction_functions[j](species_, temperature_, log_temperature_, pressure_, mixture_concentration_);
                    auto point_progress = progress_rate_functions[j](species_, temperature_, gibbs_free_energies_, point_forward_reaction);
                    accrue_source_on_fly[j](point_progress, point_source_local[i - start]);
                }}
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
        self.write_reaction_functions(file, reaction_calls, configuration)
        self.write_progress_rate_functions(file, progress_rates, is_reversible, equilibrium_constants, configuration)
        self.write_species_production_functions_on_the_fly(file, species_production_on_fly_function_texts, configuration)

        self.write_start_of_source_function_threaded(file, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
        self.write_reaction_functions_pointer_list(file, reaction_calls, configuration)
        self.write_progress_rate_functions_pointer_list(file, reaction_calls, configuration)
        self.write_accrue_source_pointer_list(file, species_production_on_fly_function_texts, configuration)
        self.write_point_loop(file, configuration)

        self.write_end_of_function(file, configuration)
        headers.append('source.h')