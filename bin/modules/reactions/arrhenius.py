def arrhenius_text(i, A, B, E, configuration):
    return_text = "{device_option}\n{scalar_function} call_forward_reaction_{i}({scalar_parameter} temperature, {scalar_parameter} log_temperature) {const_option} {{ return arrhenius({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature);}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B)

def darrhenius_text_dtemperature(i, A, B, E, configuration):
    return_text = "{device_option}\n{scalar_function} dcall_forward_reaction_{i}_dtemperature({scalar_parameter} temperature, {scalar_parameter} log_temperature) {const_option} {{ return darrhenius_dtemperature({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature);}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B)

def darrhenius_text_dlog_temperature(i, A, B, E, configuration):
    return_text = "{device_option}\n{scalar_function} dcall_forward_reaction_{i}_dlog_temperature({scalar_parameter} temperature, {scalar_parameter} log_temperature) {const_option} {{ return darrhenius_dlog_temperature({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature);}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B)

def create_reaction_functions_and_calls_arrhenius(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False):
        if verbose:
            print(f"  Arrhenius Parameters: A = {reaction.rate.pre_exponential_factor}, "
                f"b = {reaction.rate.temperature_exponent}, "
                f"Ea = {reaction.rate.activation_energy}")
        reaction_rates[reaction_index] = arrhenius_text(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, configuration)
        reaction_rates_derivatives.append(darrhenius_text_dtemperature(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, configuration))
        reaction_rates_derivatives.append(darrhenius_text_dlog_temperature(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, configuration))
        reaction_calls[reaction_index] = "call_forward_reaction_{reaction_index}(temperature, log_temperature);\n".format(**vars(configuration),reaction_index = reaction_index)

