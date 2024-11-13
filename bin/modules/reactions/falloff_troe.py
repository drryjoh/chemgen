from .reaction_utility import *
def create_reaction_functions_and_calls_troe(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False):
    reaction_rate = reaction.rate
    if verbose:
        print(f"  Arrhenius Parameters (high pressure limit): A = {reaction_rate.high_rate.pre_exponential_factor}, "
            f"b = {reaction_rate.high_rate.temperature_exponent}, "
            f"Ea = {reaction_rate.high_rate.activation_energy}")
        print(f"  Arrhenius Parameters (low pressure limit): A = {reaction_rate.low_rate.pre_exponential_factor}, "
            f"b = {reaction_rate.low_rate.temperature_exponent}, "
            f"Ea = {reaction_rate.low_rate.activation_energy}")
    
    falloff_coeffs = reaction.rate.falloff_coeffs
    if len(falloff_coeffs) >3:
        [alpha, T3, T1, T2] = falloff_coeffs
    else:
        [alpha, T3, T1] = falloff_coeffs
        T2 = 0
    reaction_rates[reaction_index] = troe_text(reaction_index,
                                               reaction_rate.low_rate.pre_exponential_factor, reaction_rate.low_rate.temperature_exponent, reaction_rate.low_rate.activation_energy,
                                               reaction_rate.high_rate.pre_exponential_factor, reaction_rate.high_rate.temperature_exponent, reaction_rate.high_rate.activation_energy,
                                               alpha, T1, T2, T3,
                                               reaction.efficiencies, species_names,
                                               configuration)
                            
    reaction_calls[reaction_index] = " call_forward_reaction_{reaction_index}(species, temperature, log_temperature, mixture_concentration);\n".format(**vars(configuration),reaction_index = reaction_index)    

def troe_text(i, A_low, B_low, E_low, A_high, B_high, E_high, alpha, T1, T2, T3, efficiencies, species_names, configuration):
    mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = ("{device_option}\n{scalar_function}\ncall_forward_reaction_{i}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) "
                  "{const_option} {{ return falloff_troe({scalar_cast}({A_low}), {scalar_cast}({B_low}), {scalar_cast}({E_low}), {scalar_cast}({A_high}), {scalar_cast}({B_high}), {scalar_cast}({E_high}), "
                  "{scalar_cast}({alpha}), {scalar_cast}({T1}), {scalar_cast}({T2}), {scalar_cast}({T3}), temperature, log_temperature, {mixture_concentration});}}")
    return return_text.format(**vars(configuration), i=i, A_low = A_low, B_low = B_low, E_low = E_low,
                                     A_high = A_high, B_high = B_high, E_high = E_high, 
                                     alpha = alpha, T1 = T1, T2 = T2, T3 = T3,
                                     mixture_concentration = mixture_concentration)

