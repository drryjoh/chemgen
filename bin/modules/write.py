import numpy as np

def write_type_defs(file, gas, configuration):
    n_species  = gas.n_species
    n_reactions  = gas.n_reactions


    file.write("""
const {index} n_species = {n_species};
const {index} n_variables = {n_species} + 1;
const {index} n_reactions = {n_reactions};
const {index} n_order_thermo = {n_thermo_order} + 1;
const {index} n_chemical_state = {n_species} + 1;
// Using alias for the array type (for example, an array of double values)
using {species} = {species_typedef};
using Reactions = {reactions_typedef};
using TemperatureMonomial = {temperature_monomial_typedef};
using TemperatureEnergyMonomial = {temperature_energy_monomial_typedef};
using TemperatureGibbsMonomial = {temperature_gibbs_monomial_typedef};
using ThermoTable = {scalar_list}<TemperatureEnergyMonomial, n_species>;
using ChemicalState = {chemical_state_tyedef};
using {jacobian} = {jacobian_typedef};
using ReactionSpecies = std::array<std::array<{scalar}, n_species>, n_reactions>;

""".format(**vars(configuration), 
n_species = int(n_species),
n_reactions = int(n_reactions), 
temperature_energy_monomial_typedef = "{temperature_monomial_typedef}".format(**vars(configuration)).replace("n_order_thermo", "n_order_thermo + 1"),
temperature_gibbs_monomial_typedef = "{temperature_monomial_typedef}".format(**vars(configuration)).replace("n_order_thermo", "n_order_thermo + 2"),
chemical_state_tyedef = "{species_typedef}".format(**vars(configuration)).replace("n_species", "n_species + 1"))
    )

    if configuration.eigen != "":
        file.write("""
using namespace Eigen;
using Triplet_ = Eigen::Triplet<{scalar}>;""".format(**vars(configuration)))

def array1d_int_text(arr):
    return ','.join(["{e}".format(e=e) for e in arr])

def write_permuted_sparsity_pattern(file, configuration):
    sp_lu_perm = configuration.sp_lu_perm
    n = sp_lu_perm.shape[0]
    text = [None]*n
    for i in range(n):
        row = sp_lu_perm[i,:]
        text[i] = "{" + array1d_int_text(row) + "}"
    sp_lu_perm_text = ",".join(text)
    file.write("\nstatic const {jacobian} sp_lu_perm = {{ {{ {sp_lu_perm_text} }} }};".format(**vars(configuration), sp_lu_perm_text=sp_lu_perm_text))

def write_permutation_indices(file, configuration):
    if configuration.eigen:
        # how to initialize static matrix - https://stackoverflow.com/questions/31549398/c-eigen-initialize-static-matrix
        # how to initialize and apply permutation matrix - https://stackoverflow.com/questions/57858014/permute-columns-of-matrix-in-eigen
        content =  "\nstatic const {indices_type} perm_indices = ({indices_type}() << {perm_text}).finished();".format(**vars(configuration), indices_type=f"Matrix<{configuration.index}, n_variables, 1>", perm_text=array1d_int_text(configuration.perm))
        content += "\nstatic const {indices_type} inv_perm_indices = ({indices_type}() << {inv_perm_text}).finished();".format(**vars(configuration), indices_type=f"Matrix<{configuration.index}, n_variables, 1>", inv_perm_text=array1d_int_text(configuration.inv_perm))
    else:
        content =  "\n{device_option} {constexpr} {scalar_list}<{index}, n_variables> perm_indices() {const_option} {{return {{{perm_text}}};}}".format(**vars(configuration), perm_text=array1d_int_text(configuration.perm))
        content += "\n{device_option} {constexpr} {scalar_list}<{index}, n_variables> inv_perm_indices() {const_option} {{return {{{inv_perm_text}}};}}".format(**vars(configuration), inv_perm_text=array1d_int_text(configuration.inv_perm))
    file.write(content)

    file.write("""
Matrix<{index}, n_variables, 1> get_perm_indices() {{ return perm_indices; }}
Matrix<{index}, n_variables, 1> get_inv_perm_indices() {{ return inv_perm_indices; }}
""".format(**vars(configuration)))

def write_lu_perm_row_col_indices(file, configuration):
    rows, cols = np.nonzero(configuration.sp_lu_perm)
    n_nonzeros = len(rows)

    content =  "\n{device_option} static const {scalar_list}<{index}, {n_nonzeros}> lu_perm_row_indices = {{{indices}}};".format(**vars(configuration), indices=array1d_int_text(rows), n_nonzeros=n_nonzeros)
    content += "\n{device_option} static const {scalar_list}<{index}, {n_nonzeros}> lu_perm_col_indices = {{{indices}}};".format(**vars(configuration), indices=array1d_int_text(cols), n_nonzeros=n_nonzeros)

    file.write(content)

    file.write("""
               
{device_option}
{index}
get_lu_perm_row_index({index} i)
{{
    return lu_perm_row_indices[i];
}}

{device_option}
{index}
get_lu_perm_col_index({index} i)
{{
    return lu_perm_col_indices[i];
}}

#if 0
std::vector<Triplet_> lu_perm_triplets;
lu_perm_triplets.reserve({n_nonzeros});

for ({index} i = 0; i < {n_nonzeros}; i++)
{{
    lu_perm_triplets.push_back(Triplet_(get_lu_perm_row_index[i], get_lu_perm_col_index[i], 0.)); // replace 0. with actual value
}}

SparseMatrix<{scalar}> lu_perm(n_variables, n_variables);
lu_perm.setFromTriplets(lu_perm_triplets.begin(), lu_perm_triplets.end());
#endif
""".format(**vars(configuration), n_nonzeros=n_nonzeros))

def write_molecular_weights(file, molecular_weights, inv_molecular_weights, configuration):
    content = "{device_option} {constexpr} {species_function} molecular_weights() {const_option} {{return {molecular_weights};}}".format(**vars(configuration), molecular_weights = molecular_weights)
    content += "\n"+"{device_option} {constexpr} {species_function} inv_molecular_weights() {const_option} {{return {inv_molecular_weights};}}".format(**vars(configuration), inv_molecular_weights = inv_molecular_weights)
    # non-constexpr wrapper for static library
    content += "\n"+"{device_option} {species_function} inv_molecular_weights_() {const_option} {{return inv_molecular_weights();}}".format(**vars(configuration))
    file.write(content)

def write_species_names(file, species_names, configuration):
    file.write("""
    // Define the species names as a fixed-size array
    #pragma once
    #include <string>

    static constexpr {scalar_list}<const char*, {n_species}> species_names_gen()
    {{
        return {{{species_list}}};
    }}

    // Return the species name for a given index
    static {string} species_name_gen({index} index)
    {{
        constexpr auto names = species_names_gen(); // Get the list of species names use auto for now
        return names[index]; // Return the name of the requested species
    }}
    // Return the species name for a given index
    {index} species_index_gen(const char* name)
    {{
        constexpr auto names = species_names_gen(); // Get the list of species names use auto for now
        for({index} i = 0; i<n_species; i++)
        {{
            if (std::strcmp(names[i], name) == 0)
            {{
                return i;
            }}
        }}
        return -1;
    }}
    """.format(**vars(configuration), species_list = ', '.join([f"\"{name}\"" for name in species_names]), n_species = len(species_names)))

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

{device_option}
{species_function} 
d{name}_dtemperature({scalar_parameter} temperature) {const_option} 
{{
    return {name}(dtemperature_monomial_dtemperature(temperature));
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

{device_option}
{species_function} 
d{name}_dtemperature({scalar_parameter} temperature) {const_option} 
{{
    return {name}(dtemperature_energy_monomial_dtemperature(temperature));
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

{device_option}
{species_function} 
d{name}_dtemperature({scalar_parameter} temperature) {const_option} 
{{
    return {name}(dtemperature_entropy_monomial_dtemperature(temperature));
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

{device_option}
{species_function} 
d{name}_dtemperature({scalar_parameter} temperature) {const_option} 
{{
    return {name}(temperature_gibbs_monomial(temperature));
}}
    """.format(**vars(configuration), thermo_fit = thermo_fit_text, name=name)
    file.write(content)

def write_gibbs_reaction_transport_fit(file, name, thermo_fit_text, configuration):
    content ="""
{device_option}
{reactions_function} 
{name}({temperature_monomial_parameter} log_temperature_monomial_sequence) {const_option} 
{{
{thermo_fit}
}}

{device_option}
{reactions_function} 
{name}({scalar_parameter} log_temperature) {const_option} 
{{
    return {name}(temperature_monomial(log_temperature));
}}

{device_option}
{reactions_function} 
d{name}_dlog_temperature({scalar_parameter} log_temperature) {const_option} 
{{
    return {name}(dtemperature_monomial_dtemperature(log_temperature)); //functionality is the same
}}
    """.format(**vars(configuration), thermo_fit = thermo_fit_text, name=name)
    file.write(content)