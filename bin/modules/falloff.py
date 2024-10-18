import sys
def create_reaction_functions_and_calls_falloff(reaction_rates
, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names):
    print(f"details for falloff reaction {reaction_index+1}")
    falloff_type = reaction.reaction_type
    if "Troe" in falloff_type:
        print("Troe")
        create_reaction_functions_and_calls_troe(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names)
        
    '''

    elif falloff_type == "SRI":
        print(f"  This is an SRI falloff reaction with parameters: {reaction.falloff.parameters}")
    
    elif falloff_type == "Lindemann":
        print(f"  This is a Lindemann falloff reaction (no additional parameters).")
    
    else:
        print(f"  Unknown falloff model: {falloff_type}")
    '''
    #sys.exit("Exiting the program")
def troe_text(i, A_low, B_low, E_low, A_high, B_high, E_high, alpha, T1, T2, T3, configuration):
    return_text = ("{device_option}\n{scalar_function}\ncall_forward_reaction_{i}({species_parameter} species, {scalar_parameter} temperature) "
                  "{const_option} {{ return falloff_troe({scalar_cast}({A_low}), {scalar_cast}({B_low}), {scalar_cast}({E_low}), {scalar_cast}({A_high}), {scalar_cast}({B_high}), {scalar_cast}({E_high}), "
                  "{scalar_cast}({alpha}), {scalar_cast}({T1}), {scalar_cast}({T2}), {scalar_cast}({T3}), temperature, mixture_concentration(species));}}")
    return return_text.format(**vars(configuration), i=i, A_low = A_low, B_low = B_low, E_low = E_low, A_high = A_high, B_high = B_high, E_high = E_high, alpha = alpha, T1 = T1, T2 = T2, T3 = T3)

def create_reaction_functions_and_calls_troe(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names):
    reaction_rate = reaction.rate
    print(f"  Arrhenius Parameters (high pressure limit): A = {reaction_rate.high_rate.pre_exponential_factor}, "
        f"b = {reaction_rate.high_rate.temperature_exponent}, "
        f"Ea = {reaction_rate.high_rate.activation_energy}")
    print(f"  Arrhenius Parameters (low pressure limit): A = {reaction_rate.low_rate.pre_exponential_factor}, "
        f"b = {reaction_rate.low_rate.temperature_exponent}, "
        f"Ea = {reaction_rate.low_rate.activation_energy}")
    
    falloff_coeffs = reaction.rate.falloff_coeffs
    print(falloff_coeffs)
    if len(falloff_coeffs) >3:
        [alpha, T3, T1, T2] = falloff_coeffs
    else:
        [alpha, T3, T1] = falloff_coeffs
        T2 = 0
    reaction_rates[reaction_index] = troe_text(reaction_index,
                                               reaction_rate.high_rate.pre_exponential_factor, reaction_rate.high_rate.temperature_exponent, reaction_rate.high_rate.activation_energy,
                                               reaction_rate.low_rate.pre_exponential_factor, reaction_rate.low_rate.temperature_exponent, reaction_rate.low_rate.activation_energy,
                                               alpha, T1, T2, T3,
                                               configuration)
    print("hello!!!")
    print(reaction_rates[reaction_index])
    reaction_calls[reaction_index] = "{scalar} forward_reaction_{reaction_index} = call_forward_reaction_{reaction_index}(species, temperature);\n".format(**vars(configuration),reaction_index = reaction_index)    

