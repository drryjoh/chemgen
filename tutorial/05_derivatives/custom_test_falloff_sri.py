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
    /*
    reference
falloff_sri(double(2.2900000000000005e+49),
            double(-11.3), 
            double(401329280.0), 
            double(5440000000000001.0), 
            double(-1.74), 
            double(361330240.0),

            double(0.138), 
            double(-670.0),
            double(0.001),
            double(1.0), 
            double(0.0), 
            temperature, 
            log_temperature, 
            mixture_concentration)
            const double& A_low,  //constant
                          const double& B_low, //constant
                          const double& E_low,  //constant
                          const double& A_high, //constant
                          const double& B_high,  //constant
                          const double& E_high, //constant
                          const double& a, //constant
                          const double& b, //constant
                          const double& c, //constant
                          const double& d, //constant
                          const double& e,// constant
                          const double& temperature, 
                          const double& mixture_concentration) 
    */

    std::ofstream outFile("falloff_sri.txt");  // Open a file to write
    std::ostream& os = outFile;  // Alias for cleaner code
    {concentration_test}
    {scalar} temperature_ =  {temperature};

    {scalar} A1_test = {scalar_cast}(2.2900000000000005e+49);
    {scalar} B1_test = {scalar_cast}(-11.3);
    {scalar} E1_test = {scalar_cast}(401329280.0);

    {scalar} A2_test = {scalar_cast}(5440000000000001);
    {scalar} B2_test = {scalar_cast}(-1.74);
    {scalar} E2_test = {scalar_cast}(361330240.0);
    {scalar} a = 0.138; {scalar} b = -670.0; {scalar} c = {scalar_cast}(0.001); {scalar} d = {scalar_cast}(1.0); {scalar} e = {scalar_cast}(0.5);

    {scalar} mixture_concentration_test = {scalar_cast}(0.5);

    {scalar} temperature_test = {scalar_cast}(1000);
    {scalar} log_temperature_test = log_gen({scalar_cast}(1000));
    {scalar} x_test = temperature_test;
    {scalar} Pr_test = {scalar_cast}(357.688);


    
//test temperature derivative
    x_test = temperature_test;
    auto my_function_0 = [&](double x) {{return falloff_sri(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, a, b, c, d, e, x, log_temperature_test, mixture_concentration_test);}};
    

    {scalar} df_dT = dfalloff_sri_dtemperature(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, a, b, c, d, e, temperature_test, log_temperature_test, mixture_concentration_test);
    std::cout << "falloff sri derivative:  " << df_dT << std::endl;
    std::cout << "falloff sri derivative Check:  ";
    
    os << "dfalloff_sri_dtemperature: ";
    os << df_dT;
    os <<", " << derivative_checker(my_function_0, x_test, 100.0, 10) << std::endl;
    
    std::cout << std::endl;

//test f-temperature derivatice
    x_test = temperature_test;
    auto my_function_1 = [&](double x) {{return f_sri(a, b, c, d, e, x, Pr_test);}};
    

    df_dT = df_sri_dtemperature(a, b, c, d, e, temperature_test, Pr_test);
    std::cout << "f sri derivative temperature:  " << df_dT << std::endl;
    std::cout << "f sri derivative temperature check:  ";
    
    os << "dfalloff_sri_dtemperature: ";
    os << df_dT;
    os <<", " << derivative_checker(my_function_1, x_test, 100.0, 10) << std::endl;
    
    std::cout << std::endl;

//test f-temperature derivatice
    x_test = Pr_test;
    auto my_function_2 = [&](double x) {{return f_sri(a, b, c, d, e, temperature_test, x);}};
    

    df_dT = df_sri_dPr(a, b, c, d, e, temperature_test, Pr_test);
    std::cout << "f sri derivative Pr:  " << df_dT << std::endl;
    std::cout << "f sri derivative Pr Check:  ";
    
    os << "dfalloff_sri_dtemperature: ";
    os << df_dT;
    os <<", " << derivative_checker(my_function_2, x_test, 30.0, 10) << std::endl;
    
    std::cout << std::endl;

// log_tempearture test
    x_test = log_temperature_test;
    auto my_function_3 = [&](double x) {{return falloff_sri(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, a, b, c, d, e, temperature_test, x, mixture_concentration_test);}};
    

    df_dT = dfalloff_sri_dlog_temperature(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, a, b, c, d, e, temperature_test, log_temperature_test, mixture_concentration_test);
    std::cout << "falloff sri derivative (log temperature):  " << df_dT << std::endl;
    std::cout << "falloff sri derivative (log temperature) Check:  ";
    
    os << "dfalloff_sri_dlog_temperature: ";
    os << df_dT;
    os <<", " << derivative_checker(my_function_3, x_test, 0.1, 10) << std::endl;
    
    std::cout << std::endl;

// log_tempearture test
    x_test = mixture_concentration_test;
    auto my_function_4 = [&](double x) {{return falloff_sri(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, a, b, c, d, e, temperature_test, log_temperature_test, x);}};
    

    df_dT = dfalloff_sri_dmixture_concentration(A1_test, B1_test, E1_test, A2_test, B2_test, E2_test, a, b, c, d, e, temperature_test, log_temperature_test, mixture_concentration_test);
    std::cout << "falloff sri derivative (mixture concentration):  " << df_dT << std::endl;
    std::cout << "falloff sri derivative (mixture concentration) Check:  ";
    
    os << "dfalloff_sri_dmixture_concentration: ";
    os << df_dT;
    os <<", " << derivative_checker(my_function_4, x_test, 0.1, 10) << std::endl;
    
    std::cout << std::endl;
    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature))