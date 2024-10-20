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
        {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
        {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);\n""".format(**vars(configuration)))

def write_start_of_source_function_threaded(file, configuration):
    file.write("""
    {device_option}
    {species_function} source_threaded({species_parameter} species, {scalar_parameter} temperature) {const_option} 
    {{
        {species} net_production_rates = {{{scalar_cast}(0)}};
        {reactions} reactions = {{}};
        {reactions} progress_rates = {{}};
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
    for species_index, species_production in enumerate(species_production_rates):
        if species_production != '':
            file.write(f"        {configuration.source_element.format(i = species_index)} = {species_production};\n") 
        else:
            file.write(f"        //source_{species_index} has no production term\n")
    file.write("\n")

def write_species_production_functions(file, species_production_rates, configuration):
    for species_index, species_production in enumerate(species_production_rates):
        if species_production != '':
            file.write("{device_option} {scalar_function} progress_rate_{species_index}({reactions_parameter} progress_rates){{ return {species_production}}};\n".format(**vars(configuration), species_production = species_production, species_index = species_index)) 
        else:
            file.write(f"        //source_{species_index} has no production term\n")
    file.write("\n")


def write_reaction_calculations(file, reaction_calls):
    for reaction_index, reaction_call in enumerate(reaction_calls):
        file.write(f"        {reaction_call}")

def write_reaction_calculations_threaded(file, reaction_calls, configuration):
    for reaction_index, reaction_call in enumerate(reaction_calls):
        call = reaction_call.split('=')[1].replace('\n','')
        file.write("{device_option}\n{scalar_function} reaction_{reaction_index}({species_parameter} species, {scalar_parameter} temperature) {const_option} {{ return {call} }}".format(**vars(configuration), reaction_index = reaction_index, call = call ))

def write_reaction_calculations_threaded(file, reaction_calls, configuration):
    for reaction_index, reaction_call in enumerate(reaction_calls):
        call = reaction_call.split('=')[1].replace('\n','')
        file.write("{device_option}\n{scalar_function} reaction_{reaction_index}({species_parameter} species, {scalar_parameter} temperature) {const_option} {{ return {call} }}".format(**vars(configuration), reaction_index = reaction_index, call = call ))

def write_progress_rates_threaded(file, progress_rates, is_reversible, equilibrium_constants, configuration):
    for i, progress_rate in enumerate(progress_rates):
        if is_reversible[i]:
            file.write("        {device_option}\n{scalar_function} production_rate_{i}({species_parameter} species, {species_parameter} gibbs_free_energies, {scalar_parameter} forward_reaction_{i}) {const_option} {{{scalar} equilibrium_constant_{i} =  {equilibrium_constant}; return {progress_rate}}}\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration), progress_rate = progress_rate.split("=")[1]))
        else:
            file.write("        {device_option}\n{scalar_function} production_rate_{i}({species_parameter} species, {species_parameter} gibbs_free_energies, {scalar_parameter} forward_reaction_{i}) {const_option} {{return {progress_rate}}}\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration), progress_rate = progress_rate.split("=")[1]))
    file.write("\n")

def write_reaction_ttb_loop(file, array, pointer_list, parameters):
    file.write(f"""
    for (int i = 0; i < n_reactions; ++i) {{
        {array}[i] = {pointer_list}[i]({parameters});
    }}

    tbb::parallel_for(0, n_reactions, [&](int i) {{
        {array}[i] = {pointer_list}[i]({parameters});
    }});
    """)

def write_species_ttb_loop(file, array, pointer_list, parameters):
    file.write(f"""
    for (int i = 0; i < n_species; ++i) {{
        {array}[i] = {pointer_list}[i]({parameters});
    }}

    tbb::parallel_for(0, n_species, [&](int i) {{
        {array}[i] = {pointer_list}[i]({parameters});
    }});
    """)

def write_function_reactions_pointer_list(file, reaction_calls, configuration):
    indentation = '        '
    reaction_call_list = ','.join([f'\n{indentation}reaction_{k}' for k in range(len(reaction_calls))])
    file.write("""
        {scalar_list}<{scalar} (*)({species_parameter}, {scalar_parameter}), n_reactions> reaction_functions = {{{reaction_call_list} }};
""".format(**vars(configuration), reaction_call_list = reaction_call_list, n_reactions = len(reaction_calls)))

def write_function_progress_rates_pointer_list(file, reaction_calls, configuration):
    indentation = '        '
    progress_rates_call_list = ','.join([f'\n{indentation}production_rate_{k}' for k in range(len(reaction_calls))])
    file.write("""
        {scalar_list}<{scalar} (*)({species_parameter}, {species_parameter}, {scalar_parameter}), n_reactions> progress_rate_functions = {{{progress_rates_call_list}}};
""".format(**vars(configuration), progress_rates_call_list = progress_rates_call_list, n_reactions = len(reaction_calls)))

def write_end_of_function(file):
    file.write("        return net_production_rates;\n    }")


def write_source_serial(file, equilibrium_constants, reaction_calls, 
                        progress_rates, is_reversible, 
                        species_production_texts, headers, configuration): 
    #
    write_equilibrium_constants(file, equilibrium_constants, configuration)
    write_start_of_source_function(file, configuration)
    
    write_reaction_calculations(file, reaction_calls)
    write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration)
    write_species_production(file, species_production_texts, configuration)
    headers.append('source.h')
    write_end_of_function(file)

def write_source_threaded(file, equilibrium_constants, reaction_calls, 
                        progress_rates, is_reversible, 
                        species_production_function_texts, headers, configuration):
    write_reaction_calculations_threaded(file, reaction_calls, configuration)
    write_progress_rates_threaded(file, progress_rates, is_reversible, equilibrium_constants, configuration)
    write_species_production_functions(file, species_production_function_texts, configuration)
    write_start_of_source_function_threaded(file, configuration)
    write_function_reactions_pointer_list(file, reaction_calls, configuration)
    write_function_progress_rates_pointer_list(file, reaction_calls, configuration)
    write_reaction_ttb_loop(file, "reactions", "reaction_functions", "species, temperature")
    write_reaction_ttb_loop(file, "progress_rates", "progress_rate_functions", "gibbs_free_energies, species, reactions[i]")
    write_species_ttb_loop(file, "progress_rates", "progress_rate_functions", "gibbs_free_energies, species, reactions[i]")
    #write_reaction_calculations_threaded(file, reaction_calls)
    #write_progress_rates_threaded(file, progress_rates, is_reversible, equilibrium_constants, configuration)
    #write_species_production_threaded(file, species_production_texts, configuration)
    #headers.append('source.h')
    #write_end_of_function_threaded(file)