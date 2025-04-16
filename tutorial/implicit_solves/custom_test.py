
import importlib.util
import cantera as ct
from pathlib import Path
import yaml
import numpy as np

def get_test_conditions(chemical_mechanism):
    config_path = Path('test_configuration.yaml')
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
    return [temperature, pressure, species_string, False]

def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 0):
    test_file = destination_folder/test_file_name
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n#include <chrono>\n")
        file.write("#include <iostream>  // For printing the result to the console\n#include <fstream>\n")
        file.write("""
        
// Overload << operator for std::array
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    for (const auto& value : arr) 
    {
        os << value << " ";
    }
    return os;
}
        
        """)
        for header in headers:
            file.write(f"#include \"{header}\"\n")

        [temperature, pressure, species_string, random] = get_test_conditions(chemical_mechanism)
        gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations

        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 

        content = """


{index}
main()
{{
    // ----- Open CSV output files -----
    std::ofstream be_file("backward_euler.txt");
    std::ofstream rk4_file("rk4.txt");
    std::ofstream sdirk2_file("sdirk2.txt");
    std::ofstream sdirk4_file("sdirk4.txt");


    {concentration_test}
    {scalar} temperature_ =  {temperature};
    {scalar} int_energy = internal_energy_volume_specific(species, temperature_);
    {chemical_state} y_init = set_chemical_state(int_energy, species);
    {chemical_state} y = y_init;
    {scalar} dt = 1e-8;
    {scalar} simple = 1;
    {scalar} t = 0;

    be_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) be_file << " " << val;
    be_file << "\\n";
    auto be_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < 4000; i++)
    {{
        y = backwards_euler(y, dt);
        t = t + dt;
        be_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) be_file << " " << val;
        be_file << "\\n";
    }}
    auto be_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> be_duration = be_end - be_start;
    std::cout << "[Backward Euler] Time elapsed: " << be_duration.count() << " seconds" << std::endl;
    
    y = y_init;
    dt = 5e-7;
    simple = 1;
    t = 0;
    sdirk2_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) sdirk2_file << " " << val;
    sdirk2_file << "\\n";
    
    auto sdirk2_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < 80; i++)
    {{
        y = sdirk2(y, dt);
        t = t + dt;
        sdirk2_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) sdirk2_file << " " << val;
        sdirk2_file << "\\n";
    }}
    auto sdirk2_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sdirk2_duration = sdirk2_end - sdirk2_start;
    std::cout << "[SDIRK2] Time elapsed: " << sdirk2_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = 1e-8;
    simple = 1;
    t = 0;
    rk4_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) rk4_file << " " << val;
    rk4_file << "\\n";
    
    auto rk4_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < 4000; i++)
    {{
        y = rk4(y, dt);
        t = t + dt;
        rk4_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) rk4_file << " " << val;
        rk4_file << "\\n";
    }}
    auto rk4_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> rk4_duration = rk4_end - rk4_start;
    std::cout << "[RK4] Time elapsed: " << rk4_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = 2e-6;
    simple = 1;
    t = 0;
    sdirk4_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) sdirk4_file << " " << val;
    sdirk4_file << "\\n";
    
    auto sdirk4_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < 20; i++)
    {{
        y = sdirk4(y, dt);
        t = t + dt;
        sdirk4_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) sdirk4_file << " " << val;
        sdirk4_file << "\\n";
    }}
    auto sdirk4_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sdirk4_duration = sdirk4_end - sdirk4_start;
    std::cout << "[SDIRK4] Time elapsed: " << sdirk4_duration.count() << " seconds" << std::endl;

    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature))
