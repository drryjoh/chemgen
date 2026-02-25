#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>  // For printing the result to the console
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
#include "reactions.h"
#include "source.h"
#include "chemical_state_functions.h"

// Overload << operator for std::array
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << "[ ";
    for (const auto& value : arr) 
    {
        os << value << " ";
    }
    os << "]";
    return os;
}

int main() {
    std::cout << "*** ChemGen ***" <<std::endl;
    Species species  = {double(0.00014403138866115228),double(0.004742897433249695),double(0.00023342455811099587),double(0.00014010841338853804),double(0.006908597020718781)};
    double temperature_ =  5292.148562675751;
    Species result = source_species(species, temperature_);


    double pressure_return = pressure(species, temperature_);
    double int_energy = internal_energy_volume_specific(species, temperature_);
    std::cout << "temperature: " << temperature_<<std::endl;
    for(int i=0; i<10; ++i)
    {
        std::cout << "temperature_ for "<< i <<" iterations: " << (temperature_ - temperature(int_energy, species, i)) / (temperature_)<<std::endl;
    }
    // Output the result
    std::cout << "Source test result:  " << result << std::endl;
    std::cout << "Cantera test result: " <<"339.75377120745264 -2.2354028342020904e-42 679.5075424149053 -679.5075424149054 6.706208502606271e-42"<<std::endl;

    
    std::cout << "ChemGen internal energy: "<< int_energy <<std::endl;
    std::cout << "Cantera internal energy: " <<"4895791.579332682"<<std::endl;

    
    std::cout << "Chemgen species cps: " << species_specific_heat_constant_pressure_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species cps: " <<"1354.467590866743 4151.631874612709 1365.099538958753 1480.8404234593265 3410.1521767149684"<<std::endl;

    
    std::cout << "Chemgen species enthalpies: " << species_enthalpy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species enthalpies: " <<"6062979.730385566 20454067.48449527 2480779.018547156 -2164353.924913196 1228374.0980919146"<<std::endl;
    
    std::cout << "Chemgen species internal energies: " << species_internal_energy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species internal energies: " <<"4687850.928770437 19640613.15864152 909862.5103686153 -3164180.582962253 -1214111.130555376"<<std::endl;

    std::cout << "Chemgen species internal entropies: " << species_entropy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species internal entropies: " <<"9198.082233324907 14072.223136979535 10037.692209388966 8094.221075366296 17029.60227605266"<<std::endl;

    
    std::cout << "Chemgen species gibbs energy: " << species_gibbs_energy_mole_specific(temperature_) <<std::endl;
    std::cout << "Cantera species gibbs energy: " <<"-1363583184.8066225 -2921953987.0029097 -1418431424.4324622 -1980412673.1844733 -1601440022.1477287"<<std::endl;
    
    std::cout << "Chemgen species gibbs reactions : " << gibbs_reaction(log_gen(temperature_)) <<std::endl;
    std::cout << "Cantera species gibbs reactions : " <<"-126.52108523696553 5.445754977874034 -5.445754977874034"<<std::endl;

    std::cout << "Pressure: " <<pressure_return <<std::endl;
    std::cout << "Temperature Monomial at 300           : " <<temperature_monomial(double(300)) <<std::endl;
    std::cout << "Temperature Energy Monomial at 300           : " <<temperature_energy_monomial(double(300)) <<std::endl;
    std::cout << "Temperature Entropy Monomial at 300           : " <<temperature_entropy_monomial(double(300)) <<std::endl;
    std::cout << "Temperature Gibbs Monomial at 300           : " <<temperature_gibbs_monomial(double(300)) <<std::endl;
    
    return 0;
}
            