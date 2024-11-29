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

    # Handle the optional 'random' key
    random = test_conditions.get('random', None)  # Default to None if 'random' is not present
    if random:
        return [0, 0, 0, random]

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
    return [temperature, pressure, species_string, False]

def get_random_TPX(gas):
    import random
    random_X = [random.random() for _ in range(gas.n_species)]
    total = sum(random_X)

    # Normalize the numbers so they sum to 1
    random_X = [x / total for x in random_X]
    return (1000 + random.random() * 1000, 10132.5 + 1013250.0 * random.random(), random_X)


def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 1):
    test_file = destination_folder/test_file_name
    check_states = True
    if check_states:
        state_files = ["bad_state_189.npy", "bad_state_623.npy", "bad_state_6565.npy", "bad_state_6565.npy"]
        mech_file = "FFCM2_model.yaml"
        n_points = len(state_files)
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <fstream> \n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")
            if "types" in header:
                        file.write(f"const int n_points = {n_points};\n")
                        file.write("using PointScalar = std::unique_ptr<{scalar_list}<{scalar}, n_points>>;\n".format(**vars(configuration)))
                        file.write("using PointReactions = std::unique_ptr<{scalar_list}<{scalar_list}<{scalar}, n_reactions>, n_points>>;\n".format(**vars(configuration)))
                        file.write("using PointSpecies = std::unique_ptr<{scalar_list}<{scalar_list}<{scalar}, n_species>, n_points>>;\n".format(**vars(configuration)))

        #[temperature, pressure, species_string, random] = get_test_conditions(chemical_mechanism)
        point_concentrations = []
        point_temperatures = []
        point_source = []
        
        for point in range(n_points):
            if check_states:
                bad_state = np.load(state_files[point])
                gas = ct.Solution(mech_file)
                temperature_b = bad_state[0]
                pressure_b = ct.gas_constant * temperature_b * np.sum(bad_state[1:])
                X_b= bad_state[1:]/np.sum(bad_state[1:])
                gas.TPX =  temperature_b, pressure_b, X_b
            else:
                gas.TPX = get_random_TPX(gas)
            point_concentrations.append(gas.concentrations)
            point_temperatures.append(gas.T)
            point_source.append(gas.net_rates_of_progress)
        concentration_test_array = []
        point_source_test_array = []
        for point in range(n_points):
            concentration_test_array.append("(*concentration_tests)[{point}] = {{{list}}};".format(list = ','.join(["{c}".format(c=c, **vars(configuration)) for c in point_concentrations[point]]), point = point))
            point_source_test_array.append("(*cantera_sources)[{point}] = {{{list}}};".format(list = ','.join(["{c}".format(c=c, **vars(configuration)) for c in point_source[point]]), point = point))
        
        temperature_tests = "(*temperatures) = {{{list}}};".format(list = ','.join([str(temp) for temp in point_temperatures]))
        concentration_tests = ' '.join(concentration_test_array)
        point_cantera_net_production_rates = ' '.join(point_source_test_array)

        content = """
{scalar_function} safe_divide({scalar_parameter} a, {scalar_parameter} b) {const_option}
{{
    if(b == 0)
    {{
        return 0;
    }}
    else
    {{
        return a/b;
    }}
}}

{scalar_function} log10_choose({scalar_parameter} a) {const_option}
{{
    if(a == 0)
    {{
        return 0;
    }}
    else
    {{
        return log10_gen(a);
    }}
}}



void l2_norm({scalar_parameter} temperature, {reactions_parameter} result, {reactions_parameter} cantera_source, std::ofstream& file) 
{{

    // Calculate L2 norm
    {scalar} l2_norm = 0.0;
    file << temperature << ", ";
    for (size_t i = 0; i < n_reactions; ++i) 
    {{
        {scalar} l2_norm_0 =  log10_choose(std::abs(result[i])) - log10_choose(std::abs(cantera_source[i]));
        //file << l2_norm_0 << ", ";

        file << result[i] <<", "<<cantera_source[i] << ", ";

        l2_norm += l2_norm_0;
    }}
    l2_norm = std::sqrt(l2_norm);

    // Write the result to the CSV file
    file << l2_norm << std::endl;
}}

{index} main() 
    {{
    PointSpecies concentration_tests = std::make_unique<{scalar_list}<{scalar_list}<{scalar}, n_species>, n_points>>();
    PointReactions cantera_sources = std::make_unique<{scalar_list}<{scalar_list}<{scalar}, n_reactions>, n_points>>();
    PointScalar temperatures = std::make_unique<{scalar_list}<{scalar}, n_points>>();

    {concentration_tests};
    {temperature_tests};
    {point_cantera_net_production_rates};

    std::ofstream file("l2_norm_results.csv");
    if (!file.is_open()) 
    {{
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1;
    }}

    file << "temperature ";
    for ({index} i = 0; i < n_reactions; i++)
    {{
        file <<", Reaction_chemgen_" << i <<", reaction_cantera_"<<i;
    }}

    file <<", l2_norm"<< std::endl;

    // Process each point
    for ({index} i = 0; i < n_points; i++)
    {{
        {reactions} result = progress_rates((*concentration_tests)[i], (*temperatures)[i]);
        l2_norm((*temperatures)[i], result, (*cantera_sources)[i], file);
    }}

    // Close the file
    file.close();

    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
                   concentration_tests = concentration_tests, 
                   temperature_tests = temperature_tests, 
                   point_cantera_net_production_rates = point_cantera_net_production_rates,
                   n_points = n_points))
