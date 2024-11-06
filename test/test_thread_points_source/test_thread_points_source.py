import importlib.util
import cantera as ct
from pathlib import Path
import yaml

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
        file.write("#include <tbb/tbb.h> // testing tbb\n")
        file.write("#include <chrono>// testing timings\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")
            if "types" in header:
                        file.write(f"const int n_points = {n_points};")
                        file.write("using PointState = {scalar_list}<{scalar_list}<{scalar}, n_species+1>, n_points>;\n".format(**vars(configuration)))
                        file.write("using ChemicalState = {scalar_list}<{scalar}, n_species+1>;\n".format(**vars(configuration)))
        [temperature, pressure, species_string] = get_test_conditions(chemical_mechanism)
        gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations
        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 
        enthalpies = gas.standard_enthalpies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        entropies = gas.standard_entropies_R * ct.gas_constant/gas.molecular_weights
        energies = gas.standard_int_energies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        gibbs = gas.standard_gibbs_RT * gas.T * ct.gas_constant
        equilibrium_constants = gas.equilibrium_constants
        gas.TPX = temperature, pressure, species_string
        print(species_string)
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
    // Call the arrhenius function with the specified parameters
    {concentration_test}
    {scalar} temperature_ =  {temperature};


    {species} result_threaded = source_threaded(species, temperature_);

    // Output the result
    std::cout << "Source test result:  " << result_threaded << std::endl;
    std::cout << "Cantera test result: " <<"{cantera_net_production_rates}"<<std::endl;
    return 0;
}}
            """

        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature, 
        cantera_net_production_rates = ' '.join([f"{npr}" for npr in gas.net_production_rates])))