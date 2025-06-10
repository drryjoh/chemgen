//Thermo-fits specific heat 

TemperatureMonomial 
temperature_monomial(const double& temperature)  
{
    TemperatureMonomial temperature_monomial_sequence_ = {};
    temperature_monomial_sequence_[0] = 1.0; // The first term is always 1

    // Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i < n_order_thermo; ++i)
    {
        temperature_monomial_sequence_[i] = temperature_monomial_sequence_[i - 1] * temperature;  // T^i = T^(i-1) * T
    }
    
    return temperature_monomial_sequence_;
}


TemperatureMonomial 
dtemperature_monomial_dtemperature(const double& temperature)  
{
    TemperatureMonomial temperature_monomial_derivative_ = {};
    TemperatureMonomial temperature_monomial_ = temperature_monomial(temperature);
    temperature_monomial_derivative_[0] = 0.0; // The derivative of T^0 (constant 1) is 0


    // Derivative of Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i < n_order_thermo; ++i) 
    {
        temperature_monomial_derivative_[i] = i * temperature_monomial_[i - 1]; // i * T^(i-1)
    }
    
    return temperature_monomial_derivative_;
}


double 
contract(const TemperatureMonomial& temperature_monomial, const TemperatureMonomial& coefficients)  
{
    double fit = double(0);
    for (int i = 0; i < n_order_thermo; ++i)
    {
        fit += temperature_monomial[i] * coefficients[i];
    }
    
    return fit;
}


double 
contract_species(const Species& a, const Species& b)  
{
    double sum = double(0);
    for (int i = 0; i < n_species; ++i)
    {
        sum += a[i] * b[i];
    }
    
    return sum;
}

//thermo fits energy

TemperatureEnergyMonomial 
temperature_energy_monomial(const double& temperature)  
{
    TemperatureEnergyMonomial temperature_energy_monomial_sequence_ = {};
    temperature_energy_monomial_sequence_[0] = 1.0; // The first term is always 1

    // Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i < n_order_thermo + 1; ++i)
    {
        temperature_energy_monomial_sequence_[i] = temperature_energy_monomial_sequence_[i - 1] * temperature;  // T^i = T^(i-1) * T
    }
    
    return temperature_energy_monomial_sequence_;
}


TemperatureEnergyMonomial 
dtemperature_energy_monomial_dtemperature(const double& temperature)  
{
    TemperatureEnergyMonomial temperature_energy_monomial_derivative_ = {};
    TemperatureEnergyMonomial temperature_energy_monomial_ = temperature_energy_monomial(temperature);
    temperature_energy_monomial_derivative_[0] = 0.0; // The derivative of T^0 (constant 1) is 0


    // Derivative of Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i < n_order_thermo + 1; ++i) 
    {
        temperature_energy_monomial_derivative_[i] = i * temperature_energy_monomial_[i - 1]; // i * T^(i-1)
    }
    
    return temperature_energy_monomial_derivative_;
}


double 
contract(const TemperatureEnergyMonomial& temperature_energy_monomial, const TemperatureEnergyMonomial& coefficients)  
{
    double fit = double(0);
    for (int i = 0; i < n_order_thermo + 1; ++i)
    {
        fit += temperature_energy_monomial[i] * coefficients[i];
    }
    
    return fit;
}

//entropy
//thermo fits energy size is same, but ln(T) is at the end

TemperatureEnergyMonomial 
temperature_entropy_monomial(const double& temperature)  
{
    TemperatureEnergyMonomial temperature_energy_monomial_sequence_ = {};
    temperature_energy_monomial_sequence_[0] = double(1); // The first term is always 1
    temperature_energy_monomial_sequence_[n_order_thermo] = std::log(temperature); 

    // Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i < n_order_thermo; ++i)
    {
        temperature_energy_monomial_sequence_[i] = temperature_energy_monomial_sequence_[i - 1] * temperature;  // T^i = T^(i-1) * T
    }
    
    return temperature_energy_monomial_sequence_;
}


TemperatureEnergyMonomial 
dtemperature_entropy_monomial_dtemperature(const double& temperature)  
{
    TemperatureEnergyMonomial temperature_energy_monomial_derivative_ = {};
    TemperatureEnergyMonomial temperature_energy_monomial_ = temperature_entropy_monomial(temperature);
    temperature_energy_monomial_derivative_[0] = 0.0; // The derivative of T^0 (constant 1) is 0
    temperature_energy_monomial_derivative_[n_order_thermo] = inv_gen(temperature); // derivative of ln(T)


    // Derivative of Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i < n_order_thermo; ++i) 
    {
        temperature_energy_monomial_derivative_[i] = i * temperature_energy_monomial_[i - 1]; // i * T^(i-1)
    }
    
    return temperature_energy_monomial_derivative_;
}
//contract is the same for energy



TemperatureGibbsMonomial 
temperature_gibbs_monomial(const double& temperature)  
{
    TemperatureGibbsMonomial temperature_gibbs_monomial_sequence_ = {};
    temperature_gibbs_monomial_sequence_[0] = double(1); // The first term is always 1
    temperature_gibbs_monomial_sequence_[n_order_thermo + 1] = multiply(temperature, std::log(temperature)); // The first term is always 1

    // Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i <= n_order_thermo; ++i)
    {
        temperature_gibbs_monomial_sequence_[i] = temperature_gibbs_monomial_sequence_[i - 1] * temperature;  // T^i = T^(i-1) * T
    }
    
    return temperature_gibbs_monomial_sequence_;
}


TemperatureGibbsMonomial 
dtemperature_gibbs_monomial_dtemperature(const double& temperature)  
{
    TemperatureGibbsMonomial temperature_gibbs_monomial_derivative_ = {};
    TemperatureGibbsMonomial temperature_gibbs_monomial_ = temperature_gibbs_monomial(temperature);
    temperature_gibbs_monomial_derivative_[0] = 0.0; // The derivative of T^0 (constant 1) is 0
    temperature_gibbs_monomial_derivative_[n_order_thermo + 1] = 
    multiply_chain(temperature, double(1), std::log(temperature), inv_gen(temperature)); 


    // Derivative of Horner's method: Start with 1 and multiply progressively by T
    for (int i = 1; i <= n_order_thermo; ++i) 
    {
        temperature_gibbs_monomial_derivative_[i] = i * temperature_gibbs_monomial_[i - 1]; // i * T^(i-1)
    }
    
    return temperature_gibbs_monomial_derivative_;
}


double 
contract(const TemperatureGibbsMonomial& temperature_gibbs_monomial, const TemperatureGibbsMonomial& coefficients)  
{
    double fit = double(0);
    for (int i = 0; i < n_order_thermo + 2; ++i)
    {
        fit += temperature_gibbs_monomial[i] * coefficients[i];
    }
    
    return fit;
}

#include "thermotransport_fits.h"
//mixture quanities

//mixture concentration

double 
mixture_concentration(const Species& species)  
{
    return sum_gen(species);
}

//Pressure

double 
pressure(const Species& species, 
         const double& temperature)  
{
    return universal_gas_constant() * temperature * sum_gen(species);
}


double 
pressure_from_massfractions_density_temperature(const Species& massfractions,
                                                const double& density,
                                                const double& temperature)  
{
    return universal_gas_constant() * temperature * contract_species(inv_molecular_weights(), massfractions) * density;
}


double 
density_from_massfractions_pressure_temperature(const Species& massfractions,
                                                const double& pressure,
                                                const double& temperature)  
{
    return pressure/(universal_gas_constant() * temperature * contract_species(inv_molecular_weights(), massfractions));
}


Species 
concentrations_from_molefractions_pressure_temperature(const Species& molefractions,
                                                       const double& pressure,
                                                       const double& temperature)  
{
    double sum_concentrations = pressure/(universal_gas_constant()  * temperature);
    return scale_gen(sum_concentrations, molefractions);
}


Species 
concentrations_from_massfractions_pressure_temperature(const Species& massfractions,
                                                       const double& pressure,
                                                       const double& temperature)  
{

    return 
    scale_gen(density_from_massfractions_pressure_temperature(massfractions, pressure, temperature),
              (massfractions * inv_molecular_weights()));
}


double 
dpressure_dtemperature(const Species& species, 
                       const double& temperature)  
{
    return universal_gas_constant() * sum_gen(species);
}


Species 
dpressure_dspecies(const Species& species, 
                   const double& temperature)  
{
    return scale_gen(universal_gas_constant(), fill_gen(temperature));
}


double 
internal_energy_volume_specific(const Species& species, 
                                const double& temperature)  
{
    return  sum_gen(species * molecular_weights() * species_internal_energy_mass_specific(temperature));
}


double 
dinternal_energy_volume_specific_dtemperature(const Species& species, 
                                              const double& temperature)  
{
    return  sum_gen(species * molecular_weights() * species_internal_energy_mass_specific(dtemperature_energy_monomial_dtemperature(temperature)));
}


double 
temperature(const double& internal_energy, const Species& species, int max_iterations)  
{
    double temperature_guess  = double(500);
    for(int i=0; i < max_iterations; ++i)
    {
        temperature_guess = 
        temperature_guess - divide((internal_energy_volume_specific(species, temperature_guess) - internal_energy),
                                   dinternal_energy_volume_specific_dtemperature(species, temperature_guess));
    }
    
    return temperature_guess;
}


double 
temperature(const double& internal_energy, const Species& species)  
{
    double temperature_guess  = double(500);
    int max_iterations = 3;
    for(int i=0; i < max_iterations; ++i)
    {
        temperature_guess = 
        temperature_guess - divide((internal_energy_volume_specific(species, temperature_guess) - internal_energy),
                                   dinternal_energy_volume_specific_dtemperature(species, temperature_guess));
    }
    
    return temperature_guess;
}