#include <cmath>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>  // For printing the result to the console
#include <fstream>

//==================
#include <vector>
//==================


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
        
#include "types_inl.h"
#include "multiply_divide.h"
#include "pow_gen.h"
#include "exp_gen.h"
#include "array_handling.h"
#include "constants.h"
#include "thermally_perfect.h"
#include "arrhenius.h"
#include "third_body.h"
#include "falloff_troe.h"
#include "falloff_lindemann.h"
#include "falloff_sri.h"
#include "pressure_dependent_arrhenius.h"
#include "direct.h"
#include "reactions.h"
#include "source.h"
#include "chemical_state_functions.h"
#include "rk4.h"
/////////////////////////////////////////////
#include "./neural_net/mlp_1.h"
/////////////////////////////////////////////
#include "linear_solvers.h"
#include "backwards_euler.h"
#include "sdirk.h"
#include "rosenbroc.h"
#include "yass.h"

auto read_scalar_or_default = [](const YAML::Node& node, const std::string& key, double default_value) 
{
    if (node[key]) return node[key].as<double>();
    std::cerr << "[Warning] " << key << " not defined. Using default: " << default_value << "\n";
    return default_value;
};

Species read_species_from_yaml(const std::string& filename, 
                               double& temperature, 
                               double& pressure, 
                               double& dt_be,
                               double& dt_sdirk2,
                               double& dt_sdirk4,
                               double& dt_ros,
                               double& dt_yass,
                               double& dt_rk4,
                               double& end_time) 
{
    YAML::Node config = YAML::LoadFile(filename);
    YAML::Node test_conditions = config["test_conditions"];
    temperature = test_conditions["temperature"].as<double>();
    pressure = test_conditions["pressure"].as<double>();
    dt_be     = read_scalar_or_default(test_conditions, "dt_be",     5e-8);
    dt_sdirk2 = read_scalar_or_default(test_conditions, "dt_sdirk2", 5e-7);
    dt_ros = read_scalar_or_default(test_conditions, "dt_ros", 5e-7);
    dt_sdirk4 = read_scalar_or_default(test_conditions, "dt_sdirk4", 2e-6);
    dt_rk4    = read_scalar_or_default(test_conditions, "dt_rk4",    1e-8);
    dt_yass    = read_scalar_or_default(test_conditions, "dt_yass",  1e-8);
    end_time  = read_scalar_or_default(test_conditions, "end_time",  1e-5);

    Species species = {}; // Zero-initialize the entire species vector

    YAML::Node species_reader;
    int molefractions  = 0;
    int massfractions  = 0;
    if (test_conditions["MoleFraction"]) 
    {
        molefractions = 1;
        species_reader = test_conditions["MoleFraction"];
        std::cout << "Using MoleFraction\n";
    }
    else if (test_conditions["MassFraction"]) 
    {
        massfractions  = 1;
        species_reader = test_conditions["MassFraction"];
        std::cout << "Using MassFraction\n";
    } 
    else
    {
        throw std::runtime_error("Error: Neither 'MoleFraction' nor 'MassFraction' is defined in test_conditions.");
    }


    for (const auto& node : species_reader) 
    {
        std::string name = node["name"].as<std::string>();
        double value = node["value"].as<double>();

        int index = species_index_gen(name.c_str());
        if (index >= 0 && index < n_species) 
        {
            species[index] = value;
        } 
        else
        {
            std::cerr << "Warning: Species \"" << name << "\" not found in species list.\n";
        }
    }

    Species concentrations  = {};
    if (massfractions == 1)
    {
        double density = density_from_massfractions_pressure_temperature(species, pressure, temperature);
        concentrations = scale_gen(density, species * inv_molecular_weights());
    }
    else if (molefractions == 1)
    {
        concentrations = concentrations_from_molefractions_pressure_temperature(species, pressure, temperature);
    }
    else
    {
        throw std::runtime_error("Neither MoleFraction nor MassFraction were defined in test_conditions.");
    }
    return concentrations;
}


int
main()
{
    // ----- Open CSV output files -----
    std::ofstream be_file("backward_euler.txt");
    std::ofstream rk4_file("rk4.txt");
    std::ofstream sdirk2_file("sdirk2.txt");
    std::ofstream sdirk4_file("sdirk4.txt");
    std::ofstream ros_file("ros.txt");
    std::ofstream yass_file("yass.txt");

    double temperature_;
    double pressure_;
    double dt_be;
    double dt_sdirk2;
    double dt_ros;
    double dt_yass;
    double dt_sdirk4;
    double dt_rk4;
    double end_time;

    Species species = 
    read_species_from_yaml("test.yaml", temperature_, pressure_, 
                           dt_be, dt_sdirk2, dt_sdirk4, dt_ros, dt_yass, dt_rk4, 
                           end_time);
    double int_energy = internal_energy_volume_specific(species, temperature_);

    ChemicalState y_init = set_chemical_state(int_energy, species);
    ChemicalState y = y_init;
    double dt = dt_be;
    int  n_run = int(end_time/dt_be);
    double t = 0;

    //######################################################################################################################################
    be_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) be_file << " " << val;
    be_file << "\n";
    std::chrono::duration<double> NN_total_time; // FOR JAY TRYING TO TIME NN CALCULATIONS
    std::chrono::duration<double> P_total_time; // FOR JAY TRYING TO TIME NN CALCULATIONS
    auto be_start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n_run; i++)
    {

        // y = backwards_euler(y, dt);
        // t = t + dt;
        // be_file << t << " " << temperature(y);
        // for (const auto& val : get_species(y)) be_file << " " << val;
        // be_file << "\n";

        //=====================================================================
        ////////////////////////////////////////////
        //// THIS IS USED FOR TRAINING PURPOSES ////
        ////////////////////////////////////////////
        bool last_step = (i == n_run - 1);
        // y = backwards_euler(y, dt, 1e-10, 10, last_step);
        y = backwards_euler(y, dt,  NN_total_time, P_total_time, 1e-12, 10);
        t = t + dt;
        be_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) be_file << " " << val;
        be_file << "\n";
        //=====================================================================

    }
    auto be_end = std::chrono::high_resolution_clock::now();
    std::cout << "Total NN Inference Time: " << NN_total_time.count() << " seconds" << std::endl; // FOR JAY TRYING TO TIME NN CALCULATIONS
    std::cout << "Total Precondition Time: " << P_total_time.count() << " seconds" << std::endl; // FOR JAY TRYING TO TIME NN CALCULATIONS
    std::chrono::duration<double> be_duration = be_end - be_start;
    std::cout << "[Backward Euler] Time elapsed: " << be_duration.count() << " seconds" << std::endl;
    auto be_adjusted_duration = be_duration - NN_total_time - P_total_time;
    std::cout << "[Backwards Euler] Adjusted Time elapsed: " << be_adjusted_duration.count() << " seconds" << std::endl;
    //######################################################################################################################################

    y = y_init;
    dt = dt_sdirk2;
    t = 0;
    n_run = int(end_time/dt_sdirk2);
    sdirk2_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) sdirk2_file << " " << val;
    sdirk2_file << "\n";
    
    auto sdirk2_start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n_run; i++)
    {
        y = sdirk2(y, dt);
        t = t + dt;
        sdirk2_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) sdirk2_file << " " << val;
        sdirk2_file << "\n";
    }
    auto sdirk2_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sdirk2_duration = sdirk2_end - sdirk2_start;
    std::cout << "[SDIRK2] Time elapsed: " << sdirk2_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_ros;
    t = 0;
    n_run = int(end_time/dt_ros);
    ros_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) ros_file << " " << val;
    ros_file << "\n";
    
    auto ros_start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n_run; i++)
    {
        y = rosenbroc(y, dt);
        t = t + dt;
        ros_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) ros_file << " " << val;
        ros_file << "\n";
    }
    auto ros_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ros_duration = ros_end - ros_start;
    std::cout << "[ROSENBROC] Time elapsed: " << ros_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_yass;
    t = 0;
    n_run = int(end_time/dt_yass);
    yass_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) yass_file << " " << val;
    yass_file << "\n";
    
    auto yass_start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n_run; i++)
    {
        y = yass(y, dt);
        t = t + dt;
        yass_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) yass_file << " " << val;
        yass_file << "\n";
    }
    auto yass_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> yass_duration = yass_end - yass_start;
    std::cout << "[YASS] Time elapsed: " << yass_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_rk4;
    n_run = int(end_time/dt_rk4);
    t = 0;
    rk4_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) rk4_file << " " << val;
    rk4_file << "\n";
    
    auto rk4_start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n_run; i++)
    {
        y = rk4(y, dt);
        t = t + dt;
        rk4_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) rk4_file << " " << val;
        rk4_file << "\n";
    }
    auto rk4_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> rk4_duration = rk4_end - rk4_start;
    std::cout << "[RK4] Time elapsed: " << rk4_duration.count() << " seconds" << std::endl;

    y = y_init;
    dt = dt_sdirk4;
    n_run = int(end_time/dt_sdirk4);
    t = 0;
    sdirk4_file << t << " " << temperature(y);
    for (const auto& val : get_species(y)) sdirk4_file << " " << val;
    sdirk4_file << "\n";
    
    auto sdirk4_start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < n_run; i++)
    {
        y = sdirk4(y, dt);
        t = t + dt;
        sdirk4_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) sdirk4_file << " " << val;
        sdirk4_file << "\n";
    }
    auto sdirk4_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sdirk4_duration = sdirk4_end - sdirk4_start;
    std::cout << "[SDIRK4] Time elapsed: " << sdirk4_duration.count() << " seconds" << std::endl;

    return 0;
}
            