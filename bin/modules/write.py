
def write_type_defs(file, gas, configuration):
    n_species  = gas.n_species
    n_reactions  = gas.n_reactions
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
using ThermoTable = {scalar_list}<TemperatureEnergyMonomial, n_species>;

""".format(**vars(configuration), 
n_species = int(n_species),
n_reactions = int(n_reactions), 
temperature_energy_monomial_typedef = "{temperature_monomial_typedef}".format(**vars(configuration)).replace("n_order_thermo", "n_order_thermo + 1"),
temperature_gibbs_monomial_typedef = "{temperature_monomial_typedef}".format(**vars(configuration)).replace("n_order_thermo", "n_order_thermo + 2"))
    )

def write_molecular_weights(file, molecular_weights, configuration):
    #{device_option} {constexpr} {scalar_function} inv_pressure_atmosphere() {const_option} {{return {scalar_cast}(1)/{scalar_cast}(101325.0);}}
    content = "{device_option} {constexpr} {species_function} molecular_weights() {const_option} {{return {molecular_weights};}}".format(**vars(configuration), molecular_weights = molecular_weights)
    file.write(content)


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
'''
def write_energy_thermo_transport_fit_frozen_species(file, name, thermo_fit_text, configuration):
    content ="""
{device_option}
{thermo_table_function} 
({thermo_table_parameter} thermo_table, ) {const_option} 
{{
    return contract(thermo_table, temperature_energy_monomial_sequence)
}}

{device_option}
{species_function} 
internal_energy_frozen_species({thermo_table_parameter} thermo_table, {temperature_energy_monomial_parameter} temperature_energy_monomial_sequence) {const_option} 
{{
    return contract(thermo_table, temperature_energy_monomial_sequence)
}}

{device_option}
{species_function} 
internal_energy_frozen_species({thermo_table_parameter} thermo_table, {scalar_parameter} temperature) {const_option} 
{{
    return internal_energy_frozen_species(temperature_energy_monomial(temperature));
}}
    """.format(**vars(configuration))
    file.write(content)
'''
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