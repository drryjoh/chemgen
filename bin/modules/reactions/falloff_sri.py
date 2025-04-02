from .reaction_utility import *
import sys
def sri_text(i, A_low, B_low, E_low, A_high, B_high, E_high, a, b, c, d, e, efficiencies, species_names, configuration):
    mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = ("{device_option}\n{scalar_function}\ncall_forward_reaction_{i}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) "
                  "{const_option} {{ return falloff_sri({scalar_cast}({A_low}), {scalar_cast}({B_low}), {scalar_cast}({E_low}), {scalar_cast}({A_high}), {scalar_cast}({B_high}), {scalar_cast}({E_high}), "
                  "{scalar_cast}({a}), {scalar_cast}({b}), {scalar_cast}({c}), {scalar_cast}({d}), {scalar_cast}({e}), temperature, log_temperature, {mixture_concentration});}}")
    return return_text.format(**vars(configuration), i=i, A_low = A_low, B_low = B_low, E_low = E_low,
                                     A_high = A_high, B_high = B_high, E_high = E_high, 
                                     a = a, b = b, c = c, d = d, e = e,
                                     mixture_concentration = mixture_concentration)

def dsri_text_dtemperature(i, A_low, B_low, E_low, A_high, B_high, E_high, a, b, c, d, e, efficiencies, species_names, configuration):
    mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = ("{device_option}\n{scalar_function}\ndcall_forward_reaction_{i}_dtemperature({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) "
                  "{const_option} {{ return dfalloff_sri_dtemperature({scalar_cast}({A_low}), {scalar_cast}({B_low}), {scalar_cast}({E_low}), {scalar_cast}({A_high}), {scalar_cast}({B_high}), {scalar_cast}({E_high}), "
                  "{scalar_cast}({a}), {scalar_cast}({b}), {scalar_cast}({c}), {scalar_cast}({d}), {scalar_cast}({e}), temperature, log_temperature, {mixture_concentration});}}")
    return return_text.format(**vars(configuration), i=i, A_low = A_low, B_low = B_low, E_low = E_low,
                                     A_high = A_high, B_high = B_high, E_high = E_high, 
                                     a = a, b = b, c = c, d = d, e = e,
                                     mixture_concentration = mixture_concentration)

def dsri_text_dlog_temperature(i, A_low, B_low, E_low, A_high, B_high, E_high, a, b, c, d, e, efficiencies, species_names, configuration):
    mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = ("{device_option}\n{scalar_function}\ndcall_forward_reaction_{i}_dlog_temperature({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) "
                  "{const_option} {{ return dfalloff_sri_dlog_temperature({scalar_cast}({A_low}), {scalar_cast}({B_low}), {scalar_cast}({E_low}), {scalar_cast}({A_high}), {scalar_cast}({B_high}), {scalar_cast}({E_high}), "
                  "{scalar_cast}({a}), {scalar_cast}({b}), {scalar_cast}({c}), {scalar_cast}({d}), {scalar_cast}({e}), temperature, log_temperature, {mixture_concentration});}}")
    return return_text.format(**vars(configuration), i=i, A_low = A_low, B_low = B_low, E_low = E_low,
                                     A_high = A_high, B_high = B_high, E_high = E_high, 
                                     a = a, b = b, c = c, d = d, e = e,
                                     mixture_concentration = mixture_concentration)

def dsri_text_dmixture_concentration(i, A_low, B_low, E_low, A_high, B_high, E_high, a, b, c, d, e, efficiencies, species_names, configuration):
    mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    dmixture_concentration_dspecies = get_mixture_concentration_derivatives(efficiencies, species_names, configuration)
    return_text = ("""{device_option}\n{species_function}\ndcall_forward_reaction_{i}_dspecies({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) {const_option} 
    {{ 
        {species} dmixture_concentration_dspecies = {{{dmixture_concentration_dspecies}}};
        return scale_gen(dfalloff_sri_dmixture_concentration({scalar_cast}({A_low}), {scalar_cast}({B_low}), {scalar_cast}({E_low}), {scalar_cast}({A_high}), {scalar_cast}({B_high}), {scalar_cast}({E_high}), {scalar_cast}({a}), {scalar_cast}({b}), {scalar_cast}({c}), {scalar_cast}({d}), {scalar_cast}({e}), temperature, log_temperature, {mixture_concentration}), dmixture_concentration_dspecies );
    }}""")
    return return_text.format(**vars(configuration), i=i, A_low = A_low, B_low = B_low, E_low = E_low,
                                     A_high = A_high, B_high = B_high, E_high = E_high, 
                                     a = a, b = b, c = c, d = d, e = e,
                                     mixture_concentration = mixture_concentration,
                                     dmixture_concentration_dspecies = dmixture_concentration_dspecies)
def create_reaction_functions_and_calls_sri(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False, temperature_jacobian = False):
    reaction_rate = reaction.rate
    if verbose:
        print(f"  Arrhenius Parameters (high pressure limit): A = {reaction_rate.high_rate.pre_exponential_factor}, "
            f"b = {reaction_rate.high_rate.temperature_exponent}, "
            f"Ea = {reaction_rate.high_rate.activation_energy}")
        print(f"  Arrhenius Parameters (low pressure limit): A = {reaction_rate.low_rate.pre_exponential_factor}, "
            f"b = {reaction_rate.low_rate.temperature_exponent}, "
            f"Ea = {reaction_rate.low_rate.activation_energy}")
    
    if reaction.third_body_name !='M':
        requires_mixture_concentration[reaction_index] = False
        third_body_multiplier = "species[{index}]".format(index = species_names.index(reaction.third_body_name))
    else:
        requires_mixture_concentration[reaction_index] = True
        third_body_multiplier = 'mixture_concentration'
    
    falloff_coeffs = reaction.rate.falloff_coeffs
    
    if len(falloff_coeffs) >3:
        [a, b, c, d, e] = falloff_coeffs
    else:
        [a, b, c, d, e] = falloff_coeffs

    reaction_rates[reaction_index] = sri_text(reaction_index,
                                               reaction_rate.low_rate.pre_exponential_factor, reaction_rate.low_rate.temperature_exponent, reaction_rate.low_rate.activation_energy,
                                               reaction_rate.high_rate.pre_exponential_factor, reaction_rate.high_rate.temperature_exponent, reaction_rate.high_rate.activation_energy,
                                               a, b, c, d, e,
                                               get_efficiencies(reaction), species_names,
                                               configuration)
    if temperature_jacobian:
        reaction_rates_derivatives.append(dsri_text_dtemperature(reaction_index,
                                                reaction_rate.low_rate.pre_exponential_factor, reaction_rate.low_rate.temperature_exponent, reaction_rate.low_rate.activation_energy,
                                                reaction_rate.high_rate.pre_exponential_factor, reaction_rate.high_rate.temperature_exponent, reaction_rate.high_rate.activation_energy,
                                                a, b, c, d, e,
                                                reaction.efficiencies, species_names,
                                                configuration))

        reaction_rates_derivatives.append(dsri_text_dlog_temperature(reaction_index,
                                                reaction_rate.low_rate.pre_exponential_factor, reaction_rate.low_rate.temperature_exponent, reaction_rate.low_rate.activation_energy,
                                                reaction_rate.high_rate.pre_exponential_factor, reaction_rate.high_rate.temperature_exponent, reaction_rate.high_rate.activation_energy,
                                                a, b, c, d, e,
                                                reaction.efficiencies, species_names,
                                                configuration))
    else:
        reaction_rates_derivatives.append(f'//dcall_forward_reaction_{reaction_index} temperature derivative unused')
    reaction_rates_derivatives.append(dsri_text_dmixture_concentration(reaction_index,
                                               reaction_rate.low_rate.pre_exponential_factor, reaction_rate.low_rate.temperature_exponent, reaction_rate.low_rate.activation_energy,
                                               reaction_rate.high_rate.pre_exponential_factor, reaction_rate.high_rate.temperature_exponent, reaction_rate.high_rate.activation_energy,
                                               a, b, c, d, e,
                                               reaction.efficiencies, species_names,
                                               configuration))
                            
    reaction_calls[reaction_index] = " call_forward_reaction_{reaction_index}(species, temperature, log_temperature, {third_body_multiplier});\n".format(**vars(configuration),reaction_index = reaction_index, third_body_multiplier = third_body_multiplier)    
                                     