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
    # Define species list
    species_list = gas.species_names

    # Define major species
    major_species = { 'H2', 'O2', 'H2O', 'AR', 'N2', 'CO', 'CO2', 'CH4', 'C2H4', 'C4H10'}
    minor_minor = {"NC3H7", "IC3H7", "C3H6", "C3H5", "CH3CCH2", "AC3H4", "PC3H4", "C3H3", "C2H5CHO", "CH3COCH3",
                   "CH3COCH2", "C2H3CHO", "C3H5OH", "NC3H7O2", "NC3H7OOH", "IC3H7O2", "IC3H7OOH", "C4H2",
                   "NC4H3", "IC4H3", "C4H4", "NC4H5", "IC4H5", "C4H5-2", "C4H6", "C4H612", "C4H6-2", "C4H7",
                   "IC4H7", "IC4H7-1", "C4H81", "C4H82", "IC4H8", "NC4H9", "SC4H9", "IC4H9", "TC4H9", "C4H10",
                   "IC4H10", "H2C4O", "CH2CHCHCHO", "CH3CHCHCO", "CH3CHCHCHO", "C3H7CHO", "IC3H7CHO",
                   "C2H5COCH3", "C2H3COCH3", "OH*", "CH*"}
    


    # Initialize the species array
    species_array = np.zeros(len(species_list))
    
    major_indices = [i for i, species in enumerate(species_list) if species in major_species]
    major_values = np.random.uniform(0.1, 0.5, len(major_indices))

    minor_indices = [i for i, species in enumerate(species_list) if (species not in major_species) and (species not in minor_minor)]
    minor_values = np.random.uniform(1e-8, 1e-6, len(minor_indices))

    minor_minor_indices = [i for i, species in enumerate(species_list) if species in minor_minor]
    minor_minor_values = np.random.uniform(1e-9, 1e-7, len(minor_minor_indices))

    # Randomly make 50% of major values zero
    #major_mask = np.random.choice([True, False], size=len(major_values), p=[0.5, 0.5])
    #major_values = major_values * major_mask

    # Randomly make 50% of minor values zero
    minor_mask = np.random.choice([True, False], size=len(minor_values), p=[0.5, 0.5])
    minor_values = minor_values * minor_mask

    minor_minor_mask = np.random.choice([True, False], size=len(minor_minor_values), p=[0.2, 0.8])
    minor_minor_values = minor_minor_values * minor_minor_mask

    for idx, value in zip(minor_indices, minor_values):
        species_array[idx] = value

    for idx, value in zip(major_indices, major_values):
        species_array[idx] = value

    for idx, value in zip(minor_minor_indices, minor_minor_values):
        species_array[idx] = value   

    species_array /= species_array.sum()
    
    return (1000 + 1000 * np.random.random(), 10132.5 + 101325.0 * 2 * np.random.random(), species_array)


def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 1):
    test_file = destination_folder/test_file_name
    with open(test_file, 'w') as file:
        file.write("#include <Kokkos_Core.hpp>\n")
        file.write("#include <cmath>\n")
        file.write("#include <fstream> \n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")
            if "types" in header:
                        file.write(f"const int n_points = {n_points};\n")
                        file.write("""using PointScalar = Kokkos::View<double*>;
using PointSpecies = Kokkos::View<double**>;
""")

        #[temperature, pressure, species_string, random] = get_test_conditions(chemical_mechanism)
        point_concentrations = []
        point_temperatures = []
        point_source = []
        for point in range(n_points):
            gas.TPX = get_random_TPX(gas) 
            point_concentrations.append(gas.concentrations)
            point_temperatures.append(gas.T)
            point_source.append(gas.net_production_rates)
        
        concentration_test_array = []
        temperature_test_array = []
        for point in range(n_points):
            temperature_test_array.append(f"h_temperatures({point}) = {point_temperatures[point]};")
            for sp, concentration in enumerate(point_concentrations[point]):
                concentration_test_array.append(f"h_concentration_tests({point},{sp}) = {concentration};")
            concentration_test_array[-1] = concentration_test_array[-1]+'\n       '
        
        temperature_tests = ' '.join(temperature_test_array)
        concentration_tests = ' '.join(concentration_test_array)


        content = """
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

{scalar_function} safe_jump({scalar_parameter} a, {scalar_parameter} b) {const_option}
{{
    if(std::abs(a) <= 1e-8 && std::abs(b))
    {{
        return 0;
    }}
    else
    {{
        return a - b;
    }}
}}

void l2_norm({scalar_parameter} temperature, {species_parameter} result, {species_parameter} cantera_source, {species_parameter} concentrations,  std::ofstream& file) 
{{

    // Calculate L2 norm
    file << temperature;
    {scalar} l2_norm = {scalar_cast}(0.0);
    for (size_t i = 0; i < n_species; ++i) 
    {{
        l2_norm += ({scalar_cast}(1)/{scalar_cast}(n_species)) * pow_gen2(safe_divide(safe_jump(result[i], cantera_source[i]), cantera_source[i]));
    }}
    file<< ", " << std::sqrt(l2_norm);
    file<< std::endl;
  
}}
void write_states({species_parameter} concentrations, {scalar_parameter} temperature,  std::ofstream& states) 
{{

    // Calculate L2 norm
    states << temperature;
    {scalar} l2_norm = 0.0;
    for (size_t i = 0; i < n_species; ++i) 
    {{
        states << ", " << concentrations[i];
    }}
    states << std::endl;
}}

{index} main(int argc, char* argv[]) {{
    std::cout << "*** ChemGen ***" <<std::endl;
        Kokkos::initialize(argc, argv);
        {{   
        PointSpecies concentration_tests_device("concentration_tests", n_points, n_species);
        PointSpecies cantera_sources_device("cantera_sources", n_points, n_species);
        PointScalar temperatures_device("temperatures", n_points);

        // Create mirror views for initialization
        auto h_concentration_tests = Kokkos::create_mirror_view(concentration_tests_device);
        auto h_temperatures = Kokkos::create_mirror_view(temperatures_device);

        {concentration_tests}
        {temperature_tests}

    Kokkos::deep_copy(concentration_tests_device, h_concentration_tests);
    Kokkos::deep_copy(temperatures_device, h_temperatures);
    Kokkos::Timer timer;

    Kokkos::parallel_for("ProcessPoints", n_points, KOKKOS_LAMBDA(int i) 
    {{
        // Extract a subview (row) of species concentrations
        auto conc_subview = Kokkos::subview(concentration_tests_device, i, Kokkos::ALL());

        // Copy subview data to a local array
        Species conc_array;  
        for (int j = 0; j < n_species; j++)
        {{
            conc_array[j] = conc_subview(j);
        }}

        // Call source() with the extracted data
        Species result = source(conc_array, temperatures_device(i));

        // Add other computations if needed...
    }});
    Kokkos::fence();

    // Stop timing and print duration
    double elapsed_time = timer.seconds();
    std::cout << "Total execution time: " << elapsed_time << " seconds" << std::endl;

    }}
    Kokkos::finalize();
    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
                   concentration_tests = concentration_tests, 
                   temperature_tests = temperature_tests, 
                   n_points = n_points))
