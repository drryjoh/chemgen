from .reaction_utility import *

def third_body_text(i, A, B, E, efficiencies, species_names, requires_mixture_concentration, configuration):
    mixture_concentration = "mixture_concentration"
    if requires_mixture_concentration:
        mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = "{device_option}\n{scalar_function} call_forward_reaction_{i}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) {const_option} {{ return third_body({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature, {mixture_concentration});}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B, mixture_concentration = mixture_concentration)

def dthird_body_dtemperature_text(i, A, B, E, configuration):
    return f"auto dreaction_dtemperature_{i}(const double& temperature) const {{ return dthird_body_dtemperature({A}, {B}, {E}, temperature)}}"

def create_reaction_functions_and_calls_third_body(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False):
    if verbose:
        print(f"  Arrhenius Parameters (3-body reaction): A = {reaction.rate.pre_exponential_factor}, "
            f"b = {reaction.rate.temperature_exponent}, "
            f"Ea = {reaction.rate.activation_energy}")
        print(f"  Collision Partner Efficiencies: {efficiencies}")
    if reaction.third_body_name !='M':
        requires_mixture_concentration[reaction_index] = False
        third_body_multiplier = "species[{index}]".format(index = species_names.index(reaction.third_body_name))
    else:
        requires_mixture_concentration[reaction_index] = True
        third_body_multiplier = 'mixture_concentration'
    reaction_rates[reaction_index] = third_body_text(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, get_efficiencies(reaction), species_names, requires_mixture_concentration[reaction_index], configuration)
    reaction_calls[reaction_index] = " call_forward_reaction_{reaction_index}(species, temperature, log_temperature, {third_body_multiplier});\n".format(**vars(configuration),reaction_index = reaction_index, third_body_multiplier = third_body_multiplier)

