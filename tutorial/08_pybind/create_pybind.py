import importlib.util
import cantera as ct
from pathlib import Path
import yaml
import numpy as np

def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 0):
    test_file = destination_folder/test_file_name
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n#include <iomanip>\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")

        content = """

struct DataPoint {{
    {scalar} temperature;
    {scalar} pressure;
    {species} massfractions;
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
        dp.temperature = std::stod(value);
        
        // Read temperature
        std::getline(ss, value, ',');
        dp.pressure = std::stod(value);
        
        // Read species concentrations into a fixed-size array
        for (size_t i = 0; i < n_species; ++i) {{
            if (!std::getline(ss, value, ',')) {{
                std::cerr << "Error: Insufficient species data in row!" << std::endl;
                return {{}};
            }}
            dp.massfractions[i] = std::stod(value);
        }}

        data.push_back(dp);
    }}
    file.close();
    return data;
}}


{index} main() {{
    std::cout << "*** ChemGen ***" <<std::endl<<std::endl;

    std::string filename = "fastest.csv"; 
    std::vector<DataPoint> dataset = readDataFromFile(filename);

    if (dataset.empty()) {{
        std::cerr << "No data loaded!" << std::endl;
        return 1;
    }}

    // Example: Process the first data point
    for (size_t i = 0; i < dataset.size(); ++i)
    {{
        const DataPoint& dp = dataset[i];
        
        {species} massfractions = dp.massfractions;
        double temperature_ = dp.temperature;
        double pressure_ = dp.pressure;

        {species} concentrations = concentrations_from_massfractions_pressure_temperature(massfractions, pressure_, temperature_);

        {jacobian} dSdy = source_jacobian(concentrations, temperature_);
        for ({index} j = 0; j<n_reactions; j++)
        {{
            std::ostringstream filename_stream;
            filename_stream << "data_fast/R_" << i <<"_"<<j<< ".txt";  // or use dp.time if preferred
            std::string jacobian_filename = filename_stream.str();
            std::ofstream jacobian_file(jacobian_filename);
            std::ostream& jac = jacobian_file;  // Alias for cleaner code
            {jacobian} dSdy_R = source_jacobian_R(concentrations, temperature_, dSdy, j);
            for({index} i = 0; i <n_species; i++)
            {{
                for({index} sp = 0; sp<n_species-1; sp++)
                {{
                    jac << dSdy_R[i][sp]<<", ";
                }}
                jac << dSdy_R[i][n_species-1] << std::endl;
            }}
            std::cout << "Saved: " << jacobian_filename << std::endl;
        }}
    }}
    return 0;
}}
            """
        file.write(content.format(**vars(configuration)))