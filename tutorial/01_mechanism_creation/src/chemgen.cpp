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
    Species species  = {double(0.0),double(0.0007522590354755689),double(0.0),double(0.0015045180709511378),double(0.0),double(0.0),double(0.0),double(0.0),double(0.0),double(0.004513554212853412),double(0.0),double(0.0),double(0.0),double(0.0)};
    double temperature_ =  1800;
    Species result = source(species, temperature_);


    double pressure_return = pressure(species, temperature_);
    double int_energy = internal_energy_volume_specific(species, temperature_);
    std::cout << "temperature: " << temperature_<<std::endl;
    for(int i=0; i<10; ++i)
    {
        std::cout << "temperature_ for "<< i <<" iterations: " << (temperature_ - temperature(int_energy, species, i)) / (temperature_)<<std::endl;
    }
    // Output the result
    std::cout << "Source test result:  " << result << std::endl;
    std::cout << "Cantera test result: " <<"0.0 -0.0025190842066153346 0.0 -0.0012595421033076673 0.0 0.0025190842066153346 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"<<std::endl;

    
    std::cout << "ChemGen internal energy: "<< int_energy <<std::endl;
    std::cout << "Cantera internal energy: " <<"231991.26384204705"<<std::endl;

    
    std::cout << "Chemgen species cps: " << species_specific_heat_constant_pressure_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species cps: " <<"20621.18704899117 16616.78802469927 1301.6863166535868 1167.0799719707766 2004.0508234338902 2781.953463285981 1650.4585771578277 2127.8710008356074 5193.1609851249505 1270.558041870319 1738.5816213185194 1281.8044880331277 1357.1758540814728 2075.9084831417763"<<std::endl;

    
    std::cout << "Chemgen species enthalpies: " << species_enthalpy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species enthalpies: " <<"247236979.76114586 22924287.42497859 17546354.348052036 1614086.106703522 4952835.895921169 -9933592.095499126 2510821.7456610794 -1281946.3720538646 7799348.825509908 1746606.8628090904 62268039.86086621 -2180056.744230426 -7139104.986552954 28129761.151496165"<<std::endl;
    
    std::cout << "Chemgen species internal energies: " << species_internal_energy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species internal energies: " <<"232389725.0858722 15500660.08734177 16610918.838790463 1146368.3520727376 4072843.38032901 -10764346.062342083 2057388.0453436875 -1721942.6298499438 4060272.9162199423 1212372.811524881 61022012.659744255 -2714367.0874177096 -7479172.534388076 27249768.635904007"<<std::endl;

    std::cout << "Chemgen species internal entropies: " << species_entropy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species internal entropies: " <<"150882.52021782557 91714.63328526431 12441.629169912661 8274.87307271925 14037.717709743141 14408.506261912842 9344.433258665218 9909.851918999964 40854.761228438154 8862.383000673015 16276.6431797554 9099.375704537313 6882.076039934863 13803.098376233454"<<std::endl;

    
    std::cout << "Chemgen species gibbs energy: " << species_gibbs_energy_mole_specific(temperature_) <<std::endl;
    std::cout << "Cantera species gibbs energy: " <<"-24546369.08398769 -286598697.8168103 -77572401.94649428 -424955372.20326775 -345498157.07935154 -646178294.1554644 -472288072.9066179 -650336789.6101968 -263127938.99676064 -397957990.63080305 396003656.55118823 -519835713.6772564 -859356783.3478972 55854118.55121097"<<std::endl;
    
    std::cout << "Chemgen species gibbs reactions : " << gibbs_reaction(log_gen(temperature_)) <<std::endl;
    std::cout << "Cantera species gibbs reactions : " <<"-19.658103528322368 -19.658103528322368 34.71281695121387 -19.658103528322368 -19.658103528322368 -19.658103528322368 -19.658103528322368"<<std::endl;

    std::cout << "Pressure: " <<pressure_return <<std::endl;
    std::cout << "Temperature Monomial at 300           : " <<temperature_monomial(double(300)) <<std::endl;
    std::cout << "Temperature Energy Monomial at 300           : " <<temperature_energy_monomial(double(300)) <<std::endl;
    std::cout << "Temperature Entropy Monomial at 300           : " <<temperature_entropy_monomial(double(300)) <<std::endl;
    std::cout << "Temperature Gibbs Monomial at 300           : " <<temperature_gibbs_monomial(double(300)) <<std::endl;
    
    return 0;
}
            