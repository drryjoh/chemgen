import importlib.util
import cantera as ct
from pathlib import Path
import yaml
import numpy as np

def get_test_conditions(chemical_mechanism):
    config_path = Path('test_configuration.yaml')
    if "simple_tb" in chemical_mechanism:
        print("**simple TB**")
        test_file = 'test_tb_configuration.yaml'
    else:
        test_file = 'test_configuration.yaml'
    config_path = Path(test_file)
    if config_path.exists():
        with config_path.open('r') as file:
            configuration = yaml.safe_load(file)
    else:
        current_dir = Path(__file__).resolve().parent
        configuration_filename = current_dir.parent/ '{0}'.format(test_file)
        print("**Test file test_configuration.yaml not found using /test/test_configuration.yaml **")
        with open(configuration_filename, 'r') as file:
            configuration = yaml.safe_load(file)
    # Extract values from the parsed YAML data
    temperature = configuration['test_conditions']['temperature']
    pressure = configuration['test_conditions']['pressure']
    species_list = configuration['test_conditions']['species']

    # Print out the parsed values
    print("Test configuration:")
    print(f"Temperature: {temperature}")
    print(f"Pressure: {pressure}")
    print("Species:")

    for species in species_list:
        print(f"  - Name: {species['name']}, MoleFraction: {species['MoleFraction']}")
    species_string  = ' '.join([f"{species['name']}:{species['MoleFraction']} " for species in species_list])
    return [temperature, pressure, species_string]

def write_headers(file, headers):
    for header in headers:
        file.write(f"#include \"{header}\"\n")

def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder):
    test_file = destination_folder/test_file_name
    n_points = 1000
    n_species = gas.n_species
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n")
        file.write("#include <chrono>// testing timings\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")
            if "types" in header:
                        file.write(f"const int n_points = {n_points};")
                        file.write("using PointState = {scalar_list}<{scalar_list}<{scalar}, n_species+1>, n_points>;\n".format(**vars(configuration)))
                        file.write("using ChemicalState = {scalar_list}<{scalar}, n_species+1>;\n".format(**vars(configuration)))
                        file.write("using PointReactions = {scalar_list}<{scalar_list}<{scalar}, n_reactions>, n_points>;\n".format(**vars(configuration)))
                        file.write("using PointSpecies = {scalar_list}<{scalar_list}<{scalar}, n_species>, n_points>;\n".format(**vars(configuration)))
        #[temperature, pressure, species_string] = get_test_conditions(chemical_mechanism)
        chemical_state = []
        for i in range(n_points):
            temperature = np.random.uniform(300, 1000)  # Random temperature between 300 and 1000
            pressure = np.random.uniform(101325 / 10, 10 * 101325)  # Random pressure between 10132.5 and 1013250
            n_species = gas.n_species 
            # Generate random mole fractions that sum to 1
            random_values = np.random.rand(n_species)  # Random values between 0 and 1 for each species
            mole_fractions = random_values / np.sum(random_values)  
            gas.TPX = temperature, pressure, mole_fractions
            concentrations = gas.concentrations
            chemical_state.append("{{{temperature},{array}}}".format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]), temperature = temperature,**vars(configuration))) 
        point_state = ',\n    '.join(chemical_state)
        content = """
// Overload << operator for std::array
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {{
    os << "[ ";
    for (const auto& value : arr) 
    {{
        os << value << " ";
    }}
    os << "]";
    return os;
}}

int main() {{
    
    PointState point_state = 
    {{{{
{point_state}
    }}}};

        // Measure serial execution time
        auto start_serial = std::chrono::high_resolution_clock::now();
        auto source_serial_ = source_serial(point_state);
        auto end_serial = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> serial_time = end_serial - start_serial;
        std::cout << "Serial execution time: " << serial_time.count() << " seconds"<<std::endl;
    
    return 0;
}}
            """

        file.write(content.format(**vars(configuration), point_state = point_state))