from .reaction_utility import *
import numpy as np
import sys

def pressure_dependent_arrhenius_text(reaction_index, As, Bs, Es, pressures, species_names, configuration):
    choose_text = ''
    dchoose_text_dtemperature = ''
    dchoose_text_dpressure = ''
    indentation = '        '
    scalar_cast = '{scalar_cast}'.format(**vars(configuration))
    scalar = '{scalar}'.format(**vars(configuration))
    
    previous_pressure = pressures[0]
    rates = []
    drates_dtemperature = []
    rate_at_pressure = []
    drate_at_pressure_dtemperature = []
    unique_pressures = [pressures[0]]

    for k, pressure in enumerate(pressures):
        E = Es[k]
        A = As[k]
        B = Bs[k]
        # If the pressure changes, finalize the current list of rates
        if pressure != previous_pressure:
            rates.append(rate_at_pressure)  # Save the rates for the previous pressure
            drates_dtemperature.append(drate_at_pressure_dtemperature)
            rate_at_pressure = []  # Reset for the new pressure
            drate_at_pressure_dtemperature = []  # Reset for the new pressure
            unique_pressures.append(pressure)  # Add the new pressure to unique list
        
        rate_at_pressure.append(f"arrhenius({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature)")
        drate_at_pressure_dtemperature.append(f"darrhenius_dtemperature({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature)")

        # Update the previous pressure tracker
        previous_pressure = pressure

    # Append the last batch of rates after the loop
    if rate_at_pressure:
        rates.append(rate_at_pressure)
    if drate_at_pressure_dtemperature:
        drates_dtemperature.append(drate_at_pressure_dtemperature)

    
    rates_for_write = []
    for rate in rates:
        if len(rate)>1:
            rates_for_write.append(' + '.join(rate))
        else:
            rates_for_write.append(rate[0])
    
    drates_for_write_dtemperature = []
    for rate in drates_dtemperature:
        if len(rate)>1:
            drates_for_write_dtemperature.append(' + '.join(rate))
        else:
            drates_for_write_dtemperature.append(rate[0])


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

# Temperature derivative
    for k, pressure in enumerate(unique_pressures):
        if k ==0:
            dchoose_text_dtemperature +=f"{indentation}/**/if (log_pressure < {np.log(pressure)}) {{ return {drates_for_write_dtemperature[k]}; }}\n"
        else:
            dchoose_text_dtemperature +=(f"{indentation}else if ({np.log(unique_pressures[k-1])} <= log_pressure && log_pressure < {np.log(pressure)})"
            f"\n{indentation}{{"
            f"\n{indentation}{scalar} arrhenius_1 = {rates_for_write[k-1]};"
            f"\n{indentation}{scalar} arrhenius_2 = {rates_for_write[k]};"
            f"\n{indentation}{scalar} log_k1 = log_gen(arrhenius_1);"
            f"\n{indentation}{scalar} log_k2 = log_gen(arrhenius_2); "
            f"\n{indentation}{scalar} dlog_k1_dtemperature = log_chain(arrhenius_1,  {drates_for_write_dtemperature[k-1]});"
            f"\n{indentation}{scalar} dlog_k2_dtemperature = log_chain(arrhenius_2,  {drates_for_write_dtemperature[k]}); "
            f"\n{indentation}return dpressure_dependent_arrhenius_dtemperature(log_k1, dlog_k1_dtemperature, log_k2, dlog_k2_dtemperature, log_pressure,  {np.log(unique_pressures[k-1])}, {np.log(pressure)});\n{indentation}}}\n")
    dchoose_text_dtemperature +=f"\n{indentation}else {{ return {drates_for_write_dtemperature[-1]}; }}"

# pressure derivative
    for k, pressure in enumerate(unique_pressures):
        if k ==0:
            dchoose_text_dpressure +=f"{indentation}/**/if (log_pressure < {np.log(pressure)}) {{ return {scalar_cast}(0); }}\n"
        else:
            dchoose_text_dpressure +=(f"{indentation}else if ({np.log(unique_pressures[k-1])} <= log_pressure && log_pressure < {np.log(pressure)})"
            f"\n{indentation}{{"
            f"\n{indentation}{scalar} log_k1 = log_gen({rates_for_write[k-1]});"
            f"\n{indentation}{scalar} log_k2 = log_gen({rates_for_write[k]}); "
            f"\n{indentation}return dpressure_dependent_arrhenius_dpressure(log_k1, log_k2, log_pressure, dlog_pressure_dpressure, {np.log(unique_pressures[k-1])}, {np.log(pressure)});\n{indentation}}}\n")
    dchoose_text_dpressure +=f"\n{indentation}else {{ return {scalar_cast}(0); }}"
    return_text ="""
{device_option}
{scalar_function}
call_forward_reaction_{reaction_index}({scalar_parameter} temperature, {scalar_parameter} pressure) {const_option}
{{
        {scalar} log_pressure = log_gen(pressure);
        {scalar} rate = {scalar_cast}(0);
{choose_text}
        return rate;
}}
    """
    derivative_return_text ="""
{device_option}
{scalar_function}
dcall_forward_reaction_{reaction_index}_dtemperature({scalar_parameter} temperature, {scalar_parameter} pressure) {const_option}
{{
        {scalar} log_pressure = log_gen(pressure);
        {scalar} rate = {scalar_cast}(0);
{dchoose_text_dtemperature}
        return rate;
}}

{device_option}
{scalar_function}
dcall_forward_reaction_{reaction_index}_dpressure({scalar_parameter} temperature, {scalar_parameter} pressure) {const_option}
{{
        {scalar} inv_universal_gas_constant_temperature = inv_gen(universal_gas_constant() * temperature); 
        {scalar} log_pressure = log_gen(pressure);
        {scalar} dlog_pressure_dpressure = dlog_da(pressure);
        {scalar} rate = {scalar_cast}(0);
{dchoose_text_dpressure}
        return rate;
}}
"""
    return [return_text.format(**vars(configuration), reaction_index = reaction_index, choose_text = choose_text),
            derivative_return_text.format(**vars(configuration), reaction_index = reaction_index, dchoose_text_dpressure = dchoose_text_dpressure, dchoose_text_dtemperature = dchoose_text_dtemperature)]
                                     
def create_reaction_functions_and_calls_pressure_dependent_arrhenius(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False):
    reaction_rate = reaction.rate
    pressures = []
    As = []
    Bs = []
    Es = []

    for P, rate in reaction.rate.rates:
        if verbose:
            print(f"    At {P} Pa: A = {rate.pre_exponential_factor}, "
                f"b = {rate.temperature_exponent}, "
                f"Ea = {rate.activation_energy},"
                f"P  = {P}")
        pressures.append(P)
        As.append(rate.pre_exponential_factor)
        Bs.append(rate.temperature_exponent)
        Es.append(rate.activation_energy)

    
    [pressure_rate, pressure_rate_derivatives] = pressure_dependent_arrhenius_text(reaction_index, As, Bs, Es, pressures, species_names, configuration)
    reaction_rates[reaction_index] = pressure_rate
    reaction_rates_derivatives.append(pressure_rate_derivatives)
    reaction_calls[reaction_index] = " call_forward_reaction_{reaction_index}(temperature, pressure_);\n".format(**vars(configuration),reaction_index = reaction_index)    
