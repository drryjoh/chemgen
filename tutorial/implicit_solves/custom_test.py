
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

#include <yaml-cpp/yaml.h>


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
auto read_scalar_or_default = [](const YAML::Node& node, const std::string& key, double default_value) 
{{
    if (node[key]) return node[key].as<{scalar}>();
    std::cerr << "[Warning] " << key << " not defined. Using default: " << default_value << "\\n";
    return default_value;
}};

Species read_species_from_yaml(const std::string& filename, 
                               {scalar}& temperature, 
                               {scalar}& pressure, 
                               {scalar}& dt_be,
                               {scalar}& dt_sdirk2,
                               {scalar}& dt_sdirk4,
                               {scalar}& dt_ros,
                               {scalar}& dt_rk4,
                               {scalar}& end_time) 
{{
    YAML::Node config = YAML::LoadFile(filename);
    YAML::Node test_conditions = config["test_conditions"];
    temperature = test_conditions["temperature"].as<{scalar}>();
    pressure = test_conditions["pressure"].as<{scalar}>();
    dt_be     = read_scalar_or_default(test_conditions, "dt_be",     5e-8);
    dt_sdirk2 = read_scalar_or_default(test_conditions, "dt_sdirk2", 5e-7);
    dt_ros = read_scalar_or_default(test_conditions, "dt_ros", 5e-7);
    dt_sdirk4 = read_scalar_or_default(test_conditions, "dt_sdirk4", 2e-6);
    dt_rk4    = read_scalar_or_default(test_conditions, "dt_rk4",    1e-8);
    end_time  = read_scalar_or_default(test_conditions, "end_time",  1e-5);

    Species species = {{}}; // Zero-initialize the entire species vector

    YAML::Node species_reader;
    {index} molefractions  = 0;
    {index} massfractions  = 0;
    if (test_conditions["MoleFraction"]) 
    {{
        molefractions = 1;
        species_reader = test_conditions["MoleFraction"];
        std::cout << "Using MoleFraction\\n";
    }}
    else if (test_conditions["MassFraction"]) 
    {{
        massfractions  = 1;
        species_reader = test_conditions["MassFraction"];
        std::cout << "Using MassFraction\\n";
    }} 
    else
    {{
        throw std::runtime_error("Error: Neither 'MoleFraction' nor 'MassFraction' is defined in test_conditions.");
    }}


    for (const auto& node : species_reader) 
    {{
        std::string name = node["name"].as<std::string>();
        {scalar} value = node["value"].as<{scalar}>();

        {index} index = species_index_gen(name.c_str());
        if (index >= 0 && index < n_species) 
        {{
            species[index] = value;
        }} 
        else
        {{
            std::cerr << "Warning: Species \\"" << name << "\\" not found in species list.\\n";
        }}
    }}

    {species} concentrations  = {{}};
    if (massfractions == 1)
    {{
        {scalar} density = density_from_massfractions_pressure_temperature(species, pressure, temperature);
        concentrations = scale_gen(density, species * inv_molecular_weights());
    }}
    else if (molefractions == 1)
    {{
        concentrations = concentrations_from_molefractions_pressure_temperature(species, pressure, temperature);
    }}
    else
    {{
        throw std::runtime_error("Neither MoleFraction nor MassFraction were defined in test_conditions.");
    }}
    return concentrations;
}}


{index}
main()
{{
    // ----- Open CSV output files -----
    std::ofstream be_file("backward_euler.txt");
    std::ofstream rk4_file("rk4.txt");
    std::ofstream sdirk2_file("sdirk2.txt");
    std::ofstream sdirk4_file("sdirk4.txt");
    std::ofstream ros_file("ros.txt");

    {scalar} temperature_;
    {scalar} pressure_;
    {scalar} dt_be;
    {scalar} dt_sdirk2;
    {scalar} dt_ros;
    {scalar} dt_sdirk4;
    {scalar} dt_rk4;
    {scalar} end_time;

    Species species = 
    read_species_from_yaml("test.yaml", temperature_, pressure_, 
                           dt_be, dt_sdirk2, dt_sdirk4, dt_ros, dt_rk4, 
                           end_time);
    {scalar} int_energy = internal_energy_volume_specific(species, temperature_);

    {chemical_state} y_init = set_chemical_state(int_energy, species);
    {chemical_state} y = y_init;
    {scalar} dt = dt_be;
    {index}  n_run = {index}(end_time/dt_be);
    {scalar} t = 0;

    be_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) be_file << " " << val;
    be_file << "\\n";
    auto be_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < n_run; i++)
    {{
        y = backwards_euler(y, dt);
        t = t + dt;
        be_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) be_file << " " << val;
        be_file << "\\n";
    }}
    auto be_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<{scalar}> be_duration = be_end - be_start;
    std::cout << "[Backward Euler] Time elapsed: " << be_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_sdirk2;
    t = 0;
    n_run = {index}(end_time/dt_sdirk2);
    sdirk2_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) sdirk2_file << " " << val;
    sdirk2_file << "\\n";
    
    auto sdirk2_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < n_run; i++)
    {{
        y = sdirk2(y, dt);
        t = t + dt;
        sdirk2_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) sdirk2_file << " " << val;
        sdirk2_file << "\\n";
    }}
    auto sdirk2_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<{scalar}> sdirk2_duration = sdirk2_end - sdirk2_start;
    std::cout << "[SDIRK2] Time elapsed: " << sdirk2_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_ros;
    t = 0;
    n_run = {index}(end_time/dt_ros);
    ros_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) ros_file << " " << val;
    ros_file << "\\n";
    
    auto ros_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < n_run; i++)
    {{
        y = rosenbroc(y, dt);
        t = t + dt;
        ros_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) ros_file << " " << val;
        ros_file << "\\n";
    }}
    auto ros_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<{scalar}> ros_duration = ros_end - ros_start;
    std::cout << "[ROSENBROC] Time elapsed: " << ros_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_rk4;
    n_run = {index}(end_time/dt_rk4);
    t = 0;
    rk4_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) rk4_file << " " << val;
    rk4_file << "\\n";
    
    auto rk4_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < n_run; i++)
    {{
        y = rk4(y, dt);
        t = t + dt;
        rk4_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) rk4_file << " " << val;
        rk4_file << "\\n";
    }}
    auto rk4_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<{scalar}> rk4_duration = rk4_end - rk4_start;
    std::cout << "[RK4] Time elapsed: " << rk4_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_sdirk4;
    n_run = {index}(end_time/dt_sdirk4);
    t = 0;
    sdirk4_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) sdirk4_file << " " << val;
    sdirk4_file << "\\n";
    
    auto sdirk4_start = std::chrono::high_resolution_clock::now();
    for({index} i = 0; i < n_run; i++)
    {{
        y = sdirk4(y, dt);
        t = t + dt;
        sdirk4_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) sdirk4_file << " " << val;
        sdirk4_file << "\\n";
    }}
    auto sdirk4_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<{scalar}> sdirk4_duration = sdirk4_end - sdirk4_start;
    std::cout << "[SDIRK4] Time elapsed: " << sdirk4_duration.count() << " seconds" << std::endl;

    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature))
