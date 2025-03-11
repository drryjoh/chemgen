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
        std::cout << df_dx <<" \\n";
    }}
    std::cout<<std::endl;
    return df_dx;
}}

template <typename Func>
{temperature_monomial_function} finite_difference_monomial(Func function, {scalar_parameter} x, {scalar_parameter} dx)
{{
    return scale_gen(inv_gen(2.0 * dx),
                     function(x + dx)-function(x - dx));
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

{scalar_function} 
error_monomial({temperature_monomial_parameter} T1, 
               {temperature_monomial_parameter} T2)
{{

    {scalar} er_ = {scalar_cast}(0);
    for({index} i = 0; i < n_order_thermo; i++)
    {{
        er_ += pow_gen2(T1[i] - T2[i]);
    }}

    return pow_gen(er_,{scalar_cast}(0.5));
}}

template <typename Func>
{temperature_monomial_function} derivative_checker_monomial(Func function,  {scalar_parameter} x_test, {temperature_monomial_parameter} T1, {scalar_parameter} starting_dx, {index} n_refine)
{{
    {temperature_monomial} df_dx = {{0}};
    {scalar} err_past = {scalar_cast}(1.0);
    for({index} i = 0; i < n_refine; i++)
    {{
        {scalar} dx = starting_dx * std::pow({scalar_cast}(0.5), {scalar_cast}(i));

        df_dx = finite_difference_monomial(function, x_test, dx);
        {scalar} err_  = error_monomial(df_dx, T1);
        if(i>0)
        {{
            std::cout << (err_past/err_)/{scalar_cast}(2.0) <<std::endl;
        }}
        err_past = err_;
    }}
    std::cout<<std::endl;
    return df_dx;
}}

{index} main()
{{
    std::cout << "*** ChemGen ***" <<std::endl<<std::endl;
    std::ofstream outFile("thermo.txt");  // Open a file to write
    std::ostream& os = outFile;  // Alias for cleaner code
    {concentration_test}
    {scalar} temperature_ =  {temperature};
    {scalar} temperature_test = {scalar_cast}(1000.0);
    {scalar} x_test = temperature_test;

    auto my_function_0 = [&]({scalar} x) {{return temperature_monomial(x);}};
    {temperature_monomial} df_dT = dtemperature_monomial_dtemperature(temperature_test);
    std::cout << "temperature monomial:  " << df_dT << std::endl;
    std::cout << "temperature monomial derivative Check:\\n";
    
    os << "dtemperature_monomial_dtemperature  :" << df_dT;
    auto checker = derivative_checker_monomial(my_function_0, x_test, df_dT, 100, 10);
    os <<", " <<  checker << std::endl;

    std::cout << error_monomial(checker, df_dT);
    std::cout << std::endl;

    //species_specific_heat_constant_pressure_mass_specific
    auto my_function_1 = [&]({scalar} x) {{return species_specific_heat_constant_pressure_mass_specific(x);}};
    {species} dS_dT = dspecies_specific_heat_constant_pressure_mass_specific_dtemperature(temperature_test);
    std::cout << "Species Cp Derivative:  " << dS_dT << std::endl;
    std::cout << "temperature monomial derivative Check:\\n";
    
    os << "dspecies_specific_heat_constant_pressure_mass_specific_dtemperature  :" << dS_dT;
    auto species_checker = derivative_checker_species(my_function_1, x_test, 100, 10);
    os <<", " <<  species_checker << std::endl;

    //std::cout << error_monomial(checker, df_dT);
    std::cout << std::endl;

    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature))