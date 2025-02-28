#include "../third_party/cantera/include/cantera/thermo.h"
#include "../third_party/cantera/include/cantera/kinetics.h"
#include "../third_party/cantera/include/cantera/IdealGasMix.h"
#include <iostream>

using namespace Cantera;

int main() {
    try {
        // Create an IdealGasMix object from a chemical mechanism file (e.g., "gri30.cti" or "gri30.yaml")
        IdealGasMix gas("gri30.yaml");  // Replace "gri30.yaml" with your mechanism file

        // Set the state of the gas: temperature (in K), pressure (in Pa), and mole fractions
        gas.setState_TPX(1200.0, 101325.0, "CH4:1.0, O2:2.0, N2:7.52");

        // Get the number of species and reactions in the mechanism
        size_t n_species = gas.nSpecies();
        size_t n_reactions = gas.nReactions();

        // Create a vector to store net production rates
        std::vector<double> net_production_rates(n_species, 0.0);

        // Create a Kinetics object for the gas
        Kinetics* kinetics = gas.kinetics();

        // Compute the net production rates for each species
        kinetics->getNetProductionRates(net_production_rates.data());

        // Output the net production rates for each species
        for (size_t i = 0; i < n_species; ++i) {
            std::cout << "Species: " << gas.speciesName(i)
                      << ", Net production rate: " << net_production_rates[i] << " kmol/m^3/s" << std::endl;
        }

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }

    return 0;
}

