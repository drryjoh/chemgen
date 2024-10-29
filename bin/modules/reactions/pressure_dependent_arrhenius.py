from .reaction_utility import *
import numpy as np
import sys
def create_reaction_functions_and_calls_pressure_dependent_arrhenius(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names):
    reaction_rate = reaction.rate
    print("  PLOG Reaction with Rate Expressions:")
    pressures = []
    As = []
    Bs = []
    Es = []

    for P, rate in reaction.rate.rates:
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
    for k, [A, B, E, pressure] in enumerate(zip(As, Bs, Es, pressures)):
        if k ==0:
            choose_text +=f"{indentation}/**/ if (log_pressure < {np.log(pressure)}) {{ return arrhenius({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature); }}\n"
        else:
            choose_text +=(f"{indentation}else if ({np.log(pressures[k-1])} <= log_pressure && log_pressure < {np.log(pressure)})"
            f"\n{indentation}{{"
            f"\n{indentation}{scalar} log_k1 = log_arrhenius({np.log(As[k-1])}, {Bs[k-1]}, {Es[k-1]}, log_temperature, inv_universal_gas_constant_temperature);"
            f"\n{indentation}{scalar} log_k2 = log_arrhenius({np.log(A)}, {B}, {E}, log_temperature, inv_universal_gas_constant_temperature); "
            f"\n{indentation}return pressure_dependent_arrhenius(log_k1, log_k2, log_pressure,  {np.log(pressures[k-1])}, {np.log(pressure)});\n{indentation}}}\n")
    choose_text +=f"\n{indentation}else {{ return arrhenius({scalar_cast}({As[-1]}), {scalar_cast}({Bs[-1]}), {scalar_cast}({Es[-1]}), temperature, log_temperature); }}"

    return_text ="""
{device_option}
{scalar_function}
call_forward_reaction_{reaction_index}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} log_pressure) {const_option}
{{
        {scalar} inv_universal_gas_constant_temperature = inv(universal_gas_constant() * temperature); 
{choose_text}
}}
    """
    return return_text.format(**vars(configuration), reaction_index = reaction_index, choose_text = choose_text)
                                     
