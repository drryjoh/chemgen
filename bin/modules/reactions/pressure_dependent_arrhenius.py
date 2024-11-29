from .reaction_utility import *
import numpy as np
import sys
def create_reaction_functions_and_calls_pressure_dependent_arrhenius(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False):
    reaction_rate = reaction.rate
    pressures = []
    As = []
    Bs = []
    Es = []

    for P, rate in reaction.rate.rates:
        if verbose:
            print(f"    At {P} Pa: A = {rate.pre_exponential_factor}, "
                f"b = {rate.temperature_exponent}, "
                f"Ea = {rate.activation_energy}")
        pressures.append(P)
        As.append(rate.pre_exponential_factor)
        Bs.append(rate.temperature_exponent)
        Es.append(rate.activation_energy)


    reaction_rates[reaction_index] = pressure_dependent_arrhenius_text(reaction_index, As, Bs, Es, pressures, species_names, configuration)
    reaction_calls[reaction_index] = " call_forward_reaction_{reaction_index}(species, temperature, log_temperature, pressure_);\n".format(**vars(configuration),reaction_index = reaction_index)    

def pressure_dependent_arrhenius_text(reaction_index, As, Bs, Es, pressures, species_names, configuration):
    choose_text = ''
    indentation = '        '
    scalar_cast = '{scalar_cast}'.format(**vars(configuration))
    scalar = '{scalar}'.format(**vars(configuration))
    
    previous_pressure = pressures[0]
    rates = []
    rate_at_pressure = []
    unique_pressures = [pressures[0]]

    for k, pressure in enumerate(pressures):
        E = Es[k]
        A = As[k]
        B = Bs[k]

        # Add the current rate to the list for this pressure
        if reaction_index == 861:
            print(f"pressure: {pressure}")
            print(f"pressure_previous: {previous_pressure}")
        # If the pressure changes, finalize the current list of rates
        if pressure != previous_pressure:
            rates.append(rate_at_pressure)  # Save the rates for the previous pressure
            rate_at_pressure = []  # Reset for the new pressure
            unique_pressures.append(pressure)  # Add the new pressure to unique list
        rate_at_pressure.append(f"arrhenius({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature)")
        # Update the previous pressure tracker
        previous_pressure = pressure

    # Append the last batch of rates after the loop
    if rate_at_pressure:
        rates.append(rate_at_pressure)

    if reaction_index == 861 or reaction_index == 80 or reaction_index == 841:
        print(unique_pressures)
    
    rates_for_write = []
    for rate in rates:
        if len(rate)>1:
            rates_for_write.append(' + '.join(rate))
        else:
            rates_for_write.append(rate[0])

    for k, pressure in enumerate(unique_pressures):
        if k ==0:
            choose_text +=f"{indentation}/**/if (log_pressure < {np.log(pressure)}) {{ return {rates_for_write[k]}; }}\n"
        else:
            choose_text +=(f"{indentation}else if ({np.log(unique_pressures[k-1])} <= log_pressure && log_pressure < {np.log(pressure)})"
            f"\n{indentation}{{"
            f"\n{indentation}{scalar} log_k1 = log_gen({rates_for_write[k-1]});"
            f"\n{indentation}{scalar} log_k2 = log_gen({rates_for_write[k]}); "
            f"\n{indentation}return pressure_dependent_arrhenius(log_k1, log_k2, log_pressure,  {np.log(unique_pressures[k-1])}, {np.log(pressure)});\n{indentation}}}\n")
    choose_text +=f"\n{indentation}else {{ return {rates_for_write[-1]}; }}"

    return_text ="""
{device_option}
{scalar_function}
call_forward_reaction_{reaction_index}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} pressure) {const_option}
{{
        {scalar} inv_universal_gas_constant_temperature = inv(universal_gas_constant() * temperature); 
        {scalar} log_pressure = log_gen(pressure);
        {scalar} rate = {scalar_cast}(0);
{choose_text}
        return rate;
}}
    """
    return return_text.format(**vars(configuration), reaction_index = reaction_index, choose_text = choose_text)
                                     
