
def write_reaction_rates(file, reaction_rates):
    for reaction in reaction_rates:
        file.write(f"    {reaction}\n")

def write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration):
    for i, progress_rate in enumerate(progress_rates):
        if is_reversible[i]:
            file.write("        {scalar} equilibrium_constant_{i} = {equilibrium_constant};\n".format(i=i, equilibrium_constant = equilibrium_constants[i], **vars(configuration)))
        file.write(f"        {progress_rate}\n") 
    file.write("\n")
    
def write_species_production(file, species_production_rates, configuration = None):
    if configuration == None:
        print("Warning this may cause compilation mismatch in decorators")
        configuration = get_configuration("configuration.yaml")
    for species_index, species_production in enumerate(species_production_rates):
        if species_production != '':
            file.write(f"        {configuration.source_element.format(i = species_index)} = {species_production};\n") 
        else:
            file.write(f"        //source_{species_index} has no production term\n")
    file.write("\n")

def write_type_defs(file, gas, configuration = None):
    n_species  = gas.n_species
    n_reactions  = gas.n_reactions
    if configuration == None:
        print("Warning this may cause compilation mismatch in decorators")
        configuration = get_configuration("configuration.yaml")
    file.write("""
const {index} n_species = {n_species};
const {index} n_reactions = {n_reactions};
const {index} n_order_thermo = {n_thermo_order} + 1;
// Using alias for the array type (for example, an array of double values)
using Species = {species_typedef};
using Reactions = {reactions_typedef};
using TemperatureMonomial = {temperature_monomial_typedef};
using TemperatureEnergyMonomial = {temperature_energy_monomial_typedef};
using TemperatureGibbsMonomial = {temperature_gibbs_monomial_typedef};

""".format(**vars(configuration), 
n_species = int(n_species),
n_reactions = int(n_reactions), 
temperature_energy_monomial_typedef = "{temperature_monomial_typedef}".format(**vars(configuration)).replace("n_order_thermo", "n_order_thermo + 1"),
temperature_gibbs_monomial_typedef = "{temperature_monomial_typedef}".format(**vars(configuration)).replace("n_order_thermo", "n_order_thermo + 2"))
    )

def write_start_of_source_function(file, configuration = None):
    if configuration == None:
        print("Warning this may cause compilation mismatch in decorators")
        configuration = get_configuration("configuration.yaml")
    file.write("""
    {device_option}
    {species_function} source({species_parameter} species, {scalar_parameter} temperature) {const_option} 
    {{
        {species} net_production_rates = {{{scalar_cast}(0)}};
        {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
        {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);\n""".format(**vars(configuration)))

def write_equilibrium_constants(file, equilibrium_constants, configuration):
    equilibrium_constant_evaluations = ''
    indentation=' '*8
    print(equilibrium_constants)
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

def write_reaction_calculations(file, reaction_calls):
    for reaction_index, reaction_call in enumerate(reaction_calls):
        file.write(f"        {reaction_call}")

def write_end_of_function(file):
    file.write("        return net_production_rates;\n    }")

def write_thermo_transport_fit(file, name, thermo_fit_text, configuration):
    content ="""
{device_option}
{species_function} 
{name}({temperature_monomial_parameter} temperature_monomial_sequence) {const_option} 
{{
{thermo_fit}
}}

{device_option}
{species_function} 
{name}({scalar_parameter} temperature) {const_option} 
{{
    return {name}(temperature_monomial(temperature));
}}
    """.format(**vars(configuration), thermo_fit = thermo_fit_text, name=name)
    file.write(content)

def write_energy_thermo_transport_fit(file, name, thermo_fit_text, configuration):
    content ="""
{device_option}
{species_function} 
{name}({temperature_energy_monomial_parameter} temperature_energy_monomial_sequence) {const_option} 
{{
{thermo_fit}
}}

{device_option}
{species_function} 
{name}({scalar_parameter} temperature) {const_option} 
{{
    return {name}(temperature_energy_monomial(temperature));
}}
    """.format(**vars(configuration), thermo_fit = thermo_fit_text, name=name)
    file.write(content)

def write_entropy_thermo_transport_fit(file, name, thermo_fit_text, configuration):
    content ="""
{device_option}
{species_function} 
{name}({temperature_energy_monomial_parameter} temperature_entropy_monomial_sequence) {const_option} 
{{
{thermo_fit}
}}

{device_option}
{species_function} 
{name}({scalar_parameter} temperature) {const_option} 
{{
    return {name}(temperature_entropy_monomial(temperature));
}}
    """.format(**vars(configuration), thermo_fit = thermo_fit_text, name=name)
    file.write(content)

def write_gibbs_thermo_transport_fit(file, name, thermo_fit_text, configuration):
    content ="""
{device_option}
{species_function} 
{name}({temperature_gibbs_monomial_parameter} temperature_gibbs_monomial_sequence) {const_option} 
{{
{thermo_fit}
}}

{device_option}
{species_function} 
{name}({scalar_parameter} temperature) {const_option} 
{{
    return {name}(temperature_gibbs_monomial(temperature));
}}
    """.format(**vars(configuration), thermo_fit = thermo_fit_text, name=name)
    file.write(content)