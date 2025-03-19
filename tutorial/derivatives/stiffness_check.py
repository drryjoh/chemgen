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
        file.write("#include <sstream>  // For printing the result to the console\n#include <iomanip>\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")

        content = """

struct DataPoint {{
    double time;
    double temperature;
    Species species;
}};

// Function to read CSV file into std::array
std::vector<DataPoint> readDataFromFile(const std::string& filename) {{
    std::vector<DataPoint> data;
    std::ifstream file(filename);
    if (!file.is_open()) {{
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }}

    std::string line;
    while (std::getline(file, line)) {{
        std::stringstream ss(line);
        DataPoint dp;
        std::string value;
        
        // Read time
        std::getline(ss, value, ',');
        dp.time = std::stod(value);
        
        // Read temperature
        std::getline(ss, value, ',');
        dp.temperature = std::stod(value);
        
        // Read species concentrations into a fixed-size array
        for (size_t i = 0; i < n_species; ++i) {{
            if (!std::getline(ss, value, ',')) {{
                std::cerr << "Error: Insufficient species data in row!" << std::endl;
                return {{}};
            }}
            dp.species[i] = std::stod(value);
        }}

        data.push_back(dp);
    }}
    
    file.close();
    return data;
}}


{index} main() {{
    std::cout << "*** ChemGen ***" <<std::endl<<std::endl;
    std::ofstream outFile("simple_jacobian.txt");  // Open a file to write
    std::ostream& os = outFile;  // Alias for cleaner code

    std::string filename = "{chemical_mechanism}_data.csv"; // Replace with actual filename
    std::vector<DataPoint> dataset = readDataFromFile(filename);

    if (dataset.empty()) {{
        std::cerr << "No data loaded!" << std::endl;
        return 1;
    }}

    // Example: Process the first data point
    for (size_t i = 0; i < dataset.size(); ++i) {{
        const DataPoint& dp = dataset[i];
        
        Species species = dp.species;
        double temperature_ = dp.temperature;

        {jacobian} dSdy = source_jacobian(species, temperature_);

        std::ostringstream filename_stream;
        filename_stream << "data/jacobian_out_{chemical_mechanism}_" << i << ".txt";  // or use dp.time if preferred
        std::string jacobian_filename = filename_stream.str();
        std::ofstream jacobian_file(jacobian_filename);
        std::ostream& jac = jacobian_file;  // Alias for cleaner code
        
        {scalar} dt = {scalar_cast}(6.688963210702341e-07);
        for({index} i = 0; i <n_species; i++)
        {{
            dSdy[i][i] = dSdy[i][i]  - 1/dt;
        }}
        
        for({index} i = 0; i <n_species; i++)
        {{
            for({index} sp = 0; sp<n_species-1; sp++)
            {{
                jac << dSdy[i][sp]<<", ";
            }}
            jac << dSdy[i][n_species-1] << std::endl;
        }}
        std::cout << "Saved: " << jacobian_filename << std::endl;
    }}
    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        chemical_mechanism = chemical_mechanism))