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
    
    std::ofstream outFile("arrhenius.txt");  // Open a file to write
    std::ostream& os = outFile;  // Alias for cleaner code
    {concentration_test}
    {scalar} temperature_ =  {temperature};
    {scalar} A_test = {scalar_cast}(10000);
    {scalar} B_test = {scalar_cast}(1.5);
    {scalar} E_test = {scalar_cast}(1000);
    {scalar} temperature_test = {scalar_cast}(1000);
    {scalar} x_test = temperature_test;

//test generic expression
    x_test = 3.0;
    auto my_function_0 = [&](double x) {{return std::pow(x, {scalar_cast}(3.4));}};
    {scalar} dmy_function_0_dx = {scalar_cast}(3.4) * std::pow(x_test, {scalar_cast}(2.4));
    
    std::cout << "d/dx(x^3.4):  " << dmy_function_0_dx << std::endl;
    std::cout << "exponent in arrhenius my_function_0:  ";
    
    os << "confirm derivative checker: ";
    os << dmy_function_0_dx;
    os <<", " << derivative_checker(my_function_0, x_test, 2.0, 10) << std::endl;
    std::cout << std::endl;

//exp gen
    x_test = 3.0;
    auto my_function_0a = [&](double x) {{return exp_gen(-2.0 * x);}};

    {scalar} dmy_function_0a_dx = exp_chain(-2.0 * x_test, -2.0);
    
    os << "expgen: ";
    std::cout << "d/dx(exp(-2x)):  " << dmy_function_0a_dx << std::endl;
    std::cout << "exponent in arrhenius my_function_0a:  ";

    os << dmy_function_0a_dx;
    os <<", " << derivative_checker(my_function_0a, x_test, 2.0, 10) << std::endl;

    std::cout << std::endl;

//test exponent in arrhenius
    os << "expgen in arrhenius: ";
    x_test = 1000.0;
    auto my_function_1 = [&](double x) {{return exp_gen(divide(-E_test, universal_gas_constant() * x));}};
        
    {scalar} dexp_term_dtemperature =
        exp_chain(divide(-E_test,
                        universal_gas_constant() * temperature_test),
                ddivide_db(-E_test,
                            universal_gas_constant() * temperature_test)*universal_gas_constant());
    
    std::cout << "exponent in arrhenius Derivative:  " << dexp_term_dtemperature << std::endl;
    std::cout << "exponent in arrhenius Check:  ";
    
    os << dexp_term_dtemperature;
    os <<", " << derivative_checker(my_function_1, x_test, 10.0 , 5) << std::endl;
    
    std::cout << std::endl;

//test arrhenius
    os << "darrhenius_dtemperature(A,B,E,T): ";
    {scalar} darr_dt = darrhenius_dtemperature(A_test, B_test, E_test, temperature_test);
    auto my_function = [&](double x) {{return arrhenius(A_test, B_test, E_test, x); }};
    std::cout << "Arrhenius Derivative:  " << darr_dt << std::endl;
    std::cout << "Arrhenius Check:  ";

    os << darr_dt;
    os <<", " << derivative_checker(my_function, x_test, 10.0 , 5) << std::endl;

    std::cout << std::endl;

//test arrhenius
    os << "darrhenius_dtemperature(A,B,E,T,logT): ";
    {scalar} log_temperature_test = log_gen(temperature_test);
    {scalar} darr_dt_ = darrhenius_dtemperature(A_test, B_test, E_test, temperature_test, log_temperature_test);
    auto my_function_4 = [&](double x) {{return arrhenius(A_test, B_test, E_test, x, log_temperature_test); }};
    std::cout << "Arrhenius Derivative:  " << darr_dt_ << std::endl;
    std::cout << "Arrhenius Check:  ";

    os << darr_dt_;
    os <<", " << derivative_checker(my_function_4, x_test, 10.0 , 5) << std::endl;
  
    std::cout << std::endl;

//test arrhenius
    os << "darrhenius_dlogT(A,B,E,T,logT): ";
    x_test = log_temperature_test;
    {scalar} darr_dt_lt = darrhenius_dlog_temperature(A_test, B_test, E_test, temperature_test, log_temperature_test);
    auto my_function_5 = [&](double x) {{return arrhenius(A_test, B_test, E_test, temperature_test, x); }};
    std::cout << "Arrhenius Derivative:  " << darr_dt_lt << std::endl;
    std::cout << "Arrhenius Check:  ";

    os << darr_dt_lt;
    os <<", " << derivative_checker(my_function_5, x_test, log_gen(5.0) , 10)<< std::endl;

    std::cout << std::endl;

    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature))