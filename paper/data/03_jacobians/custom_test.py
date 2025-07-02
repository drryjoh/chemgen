import importlib.util
import cantera as ct
from pathlib import Path
import yaml
import numpy as np

def get_random_TPX(gas):
    # Define species list
    species_list = gas.species_names

    species_array =  np.random.uniform(0, 1, len(gas.species_names))
    species_array /= species_array.sum()

    return (1000 + 1500 * np.random.random(), 10132.5 + 101325.0 * 10.9 * np.random.random(), species_array)

def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 0):
    test_file = destination_folder/test_file_name
    temp_dependence = configuration.temperature_jacobian == "on"
    if temp_dependence:
        temp_dep = "1"
    else:
        temp_dep  = "0"
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n#include <iomanip>\n")
        file.write("#include <fstream>  // For printing the result to the console\n#include <iomanip>\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")

        [temperature, pressure, species_string] = get_random_TPX(gas)
        gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations

        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 

        content = """

template <typename Func>
{chemical_state_function} 
finite_difference_species_i(Func function, 
                            {chemical_state_parameter} s,
                            {chemical_state_parameter} y,
                            {scalar_parameter} dx,
                            {index} sp)
{{
    {species} x_forward = get_species(y);
    x_forward[sp]   += dx;
    {scalar} T_forward = temperature(get_energy(y), x_forward); 
    
    return scale_gen(inv_gen(dx),
                     function(x_forward, T_forward) - s);
}}

template <typename Func>
{chemical_state_function} 
finite_difference_species_i_const_T(Func function, 
                                    {chemical_state_parameter} s,
                                    {chemical_state_parameter} y,
                                    {scalar_parameter} dx,
                                    {scalar_parameter} T,
                                    {index} sp)
{{
    {species} x_forward = get_species(y);
    x_forward[sp]   += dx;
    
    return scale_gen(inv_gen(dx),
                     function(x_forward, T) - s);
}}

template <typename Func>
{chemical_state_function} 
finite_difference_temperature(Func function, 
                             {chemical_state_parameter} s,
                             {chemical_state_parameter} y,
                             {scalar_parameter} T,
                             {scalar_parameter} dT,
                             {index} sp)
{{
    {scalar} T_forward = T + dT; 
    
    return scale_gen(inv_gen(dT),
                     function(get_species(y), T_forward) - s);
}}

template <typename Func>
{chemical_state_function} 
finite_difference_T(Func function, 
                    {chemical_state_parameter} s, 
                    {chemical_state_parameter} y,
                    {scalar_parameter} T, 
                    {scalar_parameter} dT)
{{

    {scalar} T_forward = T + dT; 
    
    return scale_gen(inv_gen(dT),
                     function(get_species(y), T_forward) - s);
}}

    #include <chrono>
    #include <iostream>

    using namespace std::chrono;
{index} main() {{
    std::cout << "*** ChemGen ***" <<std::endl<<std::endl;
    {concentration_test}
    {scalar} temperature_ =  {temperature};
    {jacobian} dSdy_check = {{{scalar_cast}(0)}};
    {index} n_runs = 10000;

    {scalar} sink = {scalar_cast}(0);
    auto t0 = high_resolution_clock::now();
    for({index} it = 0; it< n_runs; it++)
    {{
        auto source_function = [&]({species} x, {scalar} T) {{return source(x, T);}};
        {chemical_state} y = set_chemical_state(internal_energy_volume_specific(species, temperature_), species);
        {chemical_state} s = source(get_species(y), temperature(y));
        for({index} sp = 0; sp < n_species; sp++)
        {{
            auto source_function = [&]({species} x, {scalar} T) {{return source(x, T);}};
#if {temp_dep}
#ifdef CHEMGEN_TEMPERATURE_EQUATION
            {chemical_state} check_dsdy =  finite_difference_species_i_const_T(source_function, s, y, 1e-6, temperature_, sp);
#else
            {chemical_state} check_dsdy =  finite_difference_species_i(source_function, s, y, 1e-6, sp);
#endif
#else
            {chemical_state} check_dsdy =  finite_difference_species_i_const_T(source_function, s, y, 1e-6, temperature_, sp);
#endif
            
            for({index} i = 0; i < n_species + 1; i++)
            {{
                dSdy_check[i][sp + 1] = check_dsdy[i];
                sink+=check_dsdy[i];
            }}
        }}
#ifdef CHEMGEN_TEMPERATURE_EQUATION
        {chemical_state} check_dsdT =  finite_difference_T(source_function, s, y, T, 1e-6);
        for({index} i = 0; i < n_species + 1; i++)
        {{
            dSdy_check[i][0] = check_dsdT[i];
            sink+=check_dsdT[i];
        }}
#endif
    }}
    auto t1 = high_resolution_clock::now();
    auto duration_fd = duration<double>(t1 - t0).count();
    std::cout << "Finite difference loop time: " << duration_fd << " ms" << std::endl;
    
    t0 = high_resolution_clock::now();
    for({index} it = 0; it< n_runs; it++)
    {{
        {jacobian} dSdy = source_jacobian(species, temperature_);
        sink+=dSdy[0][0];
    }}

    t1 = high_resolution_clock::now();
    auto duration_jac = duration<double>(t1 - t0).count();
    std::cout << "Analytical Jacobian loop time: " << duration_jac << " ms" << std::endl;
    std::cout << "Ratio: " << duration_fd/duration_jac << " ms" << std::endl;
    std::cout << sink << std::endl;
    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature,
        temp_dep = temp_dep))