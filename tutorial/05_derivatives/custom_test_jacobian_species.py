import importlib.util
import cantera as ct
from pathlib import Path
import yaml
import numpy as np

def get_test_conditions(chemical_mechanism):
    test_file = 'test_configuration.yaml'
    config_path = Path(test_file)

    if config_path.exists():
        with config_path.open('r') as file:
            configuration = yaml.safe_load(file)
    else:
        current_dir = Path(__file__).resolve().parent
        configuration_filename = current_dir.parent.parent/ 'test/{0}'.format(test_file)
        print("**Test file test_configuration.yaml not found using /test/test_configuration.yaml **")
        with open(configuration_filename, 'r') as file:
            configuration = yaml.safe_load(file)
    
    # Extract values from the parsed YAML data
    test_conditions = configuration.get('test_conditions', {})

    # Extract other required values
    temperature = test_conditions['temperature']
    pressure = test_conditions['pressure']
    species_list = test_conditions['species']

    # Print out the parsed values
    print("Test configuration:")
    print(f"Temperature: {temperature}")
    print(f"Pressure: {pressure}")
    print("Species:")

    for species in species_list:
        print(f"  - Name: {species['name']}, MoleFraction: {species['MoleFraction']}")
    species_string  = ' '.join([f"{species['name']}:{species['MoleFraction']} " for species in species_list])
    return [temperature, pressure, species_string]

def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 0):
    begin = '0'
    end = 'n_species'
    temperature_jacobian = False

    if configuration.temperature_jacobian == 'on':
        temperature_jacobian = True
        begin = '1'
        end = 'n_species + 1'
    
    test_file = destination_folder/test_file_name
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n#include <iomanip>\n")
        file.write("#include <fstream>  // For printing the result to the console\n#include <iomanip>\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")

        [temperature, pressure, species_string] = get_test_conditions(chemical_mechanism)
        gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations

        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 
        enthalpies = gas.standard_enthalpies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        entropies = gas.standard_entropies_R * ct.gas_constant/gas.molecular_weights
        energies = gas.standard_int_energies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        gibbs = gas.standard_gibbs_RT * gas.T * ct.gas_constant
        gibbs_reactions = gas.delta_standard_gibbs/gas.T/ct.gas_constant
        equilibrium_constants = gas.equilibrium_constants
        if temperature_jacobian:
            temperature_loop = """

    {scalar} x_test = temperature_;
    auto my_function_0 = [&]({scalar} x) {{return source(species, x);}};
    {species} check_dsdT =  derivative_checker_species(my_function_0, x_test, 100, 10);
    
    for({index} i = 1; i < n_species + 1; i++)
    {{
        dSdy_check[i][0] = check_dsdT[i-1];
    }}

        """.format(**vars(configuration))
        else:
            temperature_loop  = ""
        content = """

// Overload << operator for std::array
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {{
    os << std::fixed << std::setprecision(16);         // Set precision to 4 decimal places

    os << "[ ";
    for (const auto& value : arr) 
    {{
        os << value << " ";
    }}
    os << "]";
    return os;
}}

{scalar_function} safe_divide({scalar_parameter} a, {scalar_parameter} b) {const_option}
{{
    if(std::abs(b) <= 1e-10)
    {{
        return 0;
    }}
    else
    {{
        return a/b;
    }}
}}

template <typename Func>
{scalar_function} finite_difference(Func function, {scalar_parameter} x, {scalar_parameter} dx)
{{
    return (function(x + dx)-function(x - dx))/(2.0 * dx);
}}

template <typename Func>
{scalar_function} derivative_checker(Func function,  {scalar_parameter} x_test, {scalar_parameter} starting_dx, {index} n_refine)
{{
    {scalar} df_dx = {scalar_cast}(0);
    for({index} i = 0; i < n_refine; i++)
    {{
        {scalar} dx = starting_dx * std::pow({scalar_cast}(0.5), {scalar_cast}(i));

        df_dx = finite_difference(function, x_test, dx);
    }}
    return df_dx;
}}

template <typename Func>
{species_function} finite_difference_species(Func function, {scalar_parameter} x, {scalar_parameter} dx)
{{
    return scale_gen(inv_gen(2.0 * dx),
                     function(x + dx)-function(x - dx));
}}

template <typename Func>
{species_function} derivative_checker_species(Func function,  {scalar_parameter} x_test, {scalar_parameter} starting_dx, {index} n_refine)
{{
    {species} df_dx = {{0}};
    for({index} i = 0; i < n_refine; i++)
    {{
        {scalar} dx = starting_dx * std::pow({scalar_cast}(0.5), {scalar_cast}(i));
        df_dx = finite_difference_species(function, x_test, dx);
    }}
    return df_dx;
}}

template <typename Func>
{species_function} finite_difference_species_i(Func function, {species_parameter} x, {scalar_parameter} dx, {index} sp)
{{
    {species} x_forward = x;
    {species} x_backward = x;
    x_forward[sp]   = x[sp] + dx;
    x_backward[sp]  = x[sp] - dx;
    return scale_gen(inv_gen(2.0 * dx),
                     function(x_forward)-function(x_backward));
}}

template <typename Func>
{species_function} derivative_checker_species_i(Func function,  {species_parameter} x_test, {scalar_parameter} starting_dx, {index} n_refine, {index} sp)
{{
    {species} df_dx = {{0}};
    for({index} i = 0; i < n_refine; i++)
    {{
        {scalar} dx = starting_dx * std::pow({scalar_cast}(0.5), {scalar_cast}(i));
        df_dx = finite_difference_species_i(function, x_test, dx, sp);
    }}
    return df_dx;
}}

{index} main() {{
    std::cout << "*** ChemGen ***" <<std::endl<<std::endl;
    std::ofstream outFile("simple_jacobian.txt");  // Open a file to write
    std::ostream& os = outFile;  // Alias for cleaner code
    {concentration_test}
    {scalar} temperature_ =  {temperature};
    {jacobian} dSdy = source_jacobian(species, temperature_);
    {jacobian} dSdy_check = {{{scalar_cast}(0)}};

    //test all species
    for({index} sp = {begin}; sp < {end}; sp++)
    {{
        {species} species_test = species;
        auto my_function_1 = [&]({species} x) {{return source(x, temperature_);}};
        {species} check_dsdy =  derivative_checker_species_i(my_function_1, species_test, 1e-3, 10, sp);
        
        for({index} i = 0; i < n_species; i++)
        {{
            dSdy_check[i][sp] = check_dsdy[i];
        }}
    }}

{temperature_loop}

    for({index} i = 0; i < {end}; i++)
    {{
        for({index} sp = 0; sp < {end}; sp++)
        {{
            {scalar} difference  = safe_divide((dSdy[i][sp] - dSdy_check[i][sp]), dSdy[i][sp]);
            std::cout << "dSdy["<<i<<"]["<<sp<<"]= "<< dSdy[i][sp] <<", "<<dSdy_check[i][sp]<<" "<<difference<<std::endl;
        }}
    }}

    std::ofstream jacobian_file("jacobian_out.txt");  // Open a file to write
    std::ostream& jac = jacobian_file;  // Alias for cleaner code

    for({index} i = 0; i <{end}; i++)
    {{
        for({index} sp = 0; sp<{end} - 1; sp++)
        {{
            jac << dSdy[i][sp]<<", ";
        }}
        jac << dSdy[i][n_species-1] << std::endl;
    }}
    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature,
        begin = begin,
        end = end,
        temperature_loop = temperature_loop))