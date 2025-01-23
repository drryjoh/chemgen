
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

def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 0):
    test_file = destination_folder/test_file_name
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n")
        write_headers(file, headers)

        [temperature, pressure, species_string, random] = get_test_conditions(chemical_mechanism)
        gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations

        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 

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

{index} main() {{
    std::cout << "*** ChemGen ***" <<std::endl;
    {concentration_test}
    {scalar} temperature_ =  {temperature};
    {species} result = source(species, temperature_);
    {scalar} int_energy = internal_energy_volume_specific(species, temperature_);

    
    for({index} i=0; i<100; ++i)
    {{
        std::cout << "temperature_ for "<< i <<" iterations: " << (temperature_ - temperature(int_energy, species, i)) / (temperature_)<<std::endl;
    }}

    
    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature))
