def write_equilibrium_constants(file, equilibrium_constants, configuration):
    equilibrium_constant_evaluations = ''
    indentation=' '*8
    for i, equilibrium_constant_i in enumerate(equilibrium_constants):
        equilibrium_constant_evaluations += f"{indentation}equilibrium_constants[{i}] = {equilibrium_constant_i};\n"
    file.write("""
    {device_option}
    {reactions_function} equilibrium_constants({scalar_parameter} temperature) {const_option} 
    {{
        {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
        {reactions} equilibrium_constants = {{}};
        {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);
{equilibrium_constant_evaluations}
        return equilibrium_constants;
    }}
""".format(**vars(configuration),
           equilibrium_constant_evaluations = equilibrium_constant_evaluations))

def write_start_of_source_function(file, configuration):
    file.write("""
    {device_option}
    {species_function} source({species_parameter} species, {scalar_parameter} temperature) {const_option} 
    {{
        {species} net_production_rates = {{{scalar_cast}(0)}};
        {reactions} forward_reactions = {{}};
        {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
        {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);\n""".format(**vars(configuration)))

def write_reaction_rates(file, reaction_rates):
    for reaction in reaction_rates:
        file.write(f"    {reaction}\n")

def write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration):
    for i, progress_rate in enumerate(progress_rates):
        if is_reversible[i]:
            file.write("        {scalar} equilibrium_constant_{i} = {equilibrium_constant};\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration)))
        file.write(f"        {progress_rate}\n") 
    file.write("\n")
    
def write_species_production(file, species_production_rates, configuration):
    if configuration == None:
        print("Warning this may cause compilation mismatch in decorators")
        configuration = get_configuration("configuration.yaml")
    for species_index, species_production in enumerate(species_production_rates):
        if species_production != '':
            file.write(f"        {configuration.source_element.format(i = species_index)} = {species_production};\n") 
        else:
            file.write(f"        //source_{species_index} has no production term\n")
    file.write("\n")

def write_reaction_calculations(file, reaction_calls):
    for reaction_index, reaction_call in enumerate(reaction_calls):
        file.write(f"        {reaction_call}")

def write_reaction_calculations_threaded(file, reaction_calls, configuration):
    for reaction_index, reaction_call in enumerate(reaction_calls):
        call = reaction_call.split('=')[1].replace('\n','')
        file.write("{device_option}\n{scalar_function} forward_reaction_{reaction_index}({species_parameter} species, {scalar_parameter} temperature) {const_option} {{ return {call} }}".format(**vars(configuration), reaction_index = reaction_index, call = call ))

def write_reaction_ttb_loop(file):
    file.write("""
    
    for (int i = 0; i < n_reactions; ++i) {
        forward_rates[i] = forward_reaction[i](species, temperature);
    }

    tbb::parallel_for(0, n_reactions, [&](int i) {
        forward_rates[i] = forward_reaction[i](species, temperature);
    });
    """)

def write_function_pointer_list(file, reaction_calls, configuration):
    indentation = '        '
    reaction_call_list = ','.join([f'\n{indentation}forward_reaction_{k}' for k in range(len(reaction_calls))])
    file.write("""
        {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}), n_reactions> forward_reactions = {{{reaction_call_list} }};
""".format(**vars(configuration), reaction_call_list = reaction_call_list, n_reactions = len(reaction_calls)))
def write_end_of_function(file):
    file.write("        return net_production_rates;\n    }")


def write_source_serial(file, equilibrium_constants, reaction_calls, 
                        progress_rates, is_reversible, 
                        species_production_texts, headers, configuration): 
    write_reaction_calculations_threaded(file, reaction_calls, configuration)
    write_equilibrium_constants(file, equilibrium_constants, configuration)
    write_start_of_source_function(file, configuration)
    #write_function_pointer_list(file, reaction_calls, configuration)
    write_reaction_calculations(file, reaction_calls)
    write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration)
    write_species_production(file, species_production_texts, configuration)
    headers.append('source.h')
    write_end_of_function(file)

def write_source_threaded(file, equilibrium_constants, reaction_calls, 
                        progress_rates, is_reversible, 
                        species_production_texts, headers, configuration): 
    write_equilibrium_constants_threaded(file, equilibrium_constants, configuration)
    write_start_of_source_function_threaded(file, configuration)
    write_reaction_calculations_threaded(file, reaction_calls)
    write_progress_rates_threaded(file, progress_rates, is_reversible, equilibrium_constants, configuration)
    write_species_production_threaded(file, species_production_texts, configuration)
    headers.append('source.h')
    write_end_of_function_threaded(file)