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
        std::cout << df_dx <<" ";
    }}
    std::cout<<std::endl;
    return df_dx;
}}

{index} main() {{
    std::cout << "*** ChemGen ***" <<std::endl<<std::endl;
    
    std::ofstream outFile("falloff_lindemann.txt");  // Open a file to write
    std::ostream& os = outFile;  // Alias for cleaner code
    {concentration_test}
    {scalar} temperature_ =  {temperature};

    {scalar} A1_test = {scalar_cast}(1588800000000);
    {scalar} B1_test = {scalar_cast}(-2.1);
    {scalar} E1_test = {scalar_cast}(23012000);

    {scalar} A2_test = {scalar_cast}(1588800000000);
    {scalar} B2_test = {scalar_cast}(-2.1);
    {scalar} E2_test = {scalar_cast}(23012000);
    {scalar} mixture_concentration_test = {scalar_cast}(0.5);

    {scalar} temperature_test = {scalar_cast}(1000);
    {scalar} log_temperature_test = log_gen({scalar_cast}(1000));
    {scalar} x_test = temperature_test;
    /*
    referenc falloff
    falloff_lindemann(double(1588800000000.0002), 
                        double(-2.1), 
                        double(23012000.0), 
                        double(12029000.000000002),
                        double(-0.31),
                        double(29049512.0),
                        temperature, 
                        log_temperature,
                        mixture_concentration + (double(0.5)-double(1))*species[1] + (double(1.5)-double(1))*species[9]);
    */
    
//test temperature derivative
    x_test = temperature_test;
    auto my_function_0 = [&](double x) {{return falloff_lindemann(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, x, log_temperature_test, mixture_concentration_test);}};
    

    {scalar} df_dT = dfalloff_lindemann_dtemperature(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, temperature_test, log_temperature_test, mixture_concentration_test);
    std::cout << "falloff_lindemann derivative:  " << df_dT << std::endl;
    std::cout << "falloff_lindemann derivative Check:  ";
    
    os << "dfalloff_lindemann_dtemperature: ";
    os << df_dT;
    os <<", " << derivative_checker(my_function_0, x_test, 100.0, 10) << std::endl;
    
    std::cout << std::endl;

//test log_temperature derivative
    x_test = log_gen(temperature_test);
    auto my_function_1 = [&](double x) {{return falloff_lindemann(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, temperature_test, x, mixture_concentration_test);}};
    

    {scalar} df_dlogT = dfalloff_lindemann_dlog_temperature(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, temperature_test, log_temperature_test, mixture_concentration_test);
    std::cout << "falloff_lindemann derivative:  " << df_dlogT << std::endl;
    std::cout << "falloff_lindemann derivative Check:  ";
    
    os << "dfalloff_lindemann_dlog_temperature: ";
    os << df_dlogT;
    os <<", " << derivative_checker(my_function_1, x_test, log_gen(100.0), 15) << std::endl;
    
    std::cout << std::endl;
    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature))