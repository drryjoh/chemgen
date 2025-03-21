from .reaction_utility import *

def third_body_text(i, A, B, E, efficiencies, species_names, requires_mixture_concentration, configuration):
    mixture_concentration = "mixture_concentration"
    if requires_mixture_concentration:
        mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = "{device_option}\n{scalar_function} call_forward_reaction_{i}({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) {const_option} {{ return third_body({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature, {mixture_concentration});}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B, mixture_concentration = mixture_concentration)

def dthird_body_text_dtemperature(i, A, B, E, efficiencies, species_names, requires_mixture_concentration, configuration):
    mixture_concentration = "mixture_concentration"
    if requires_mixture_concentration:
        mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = "{device_option}\n{scalar_function} dcall_forward_reaction_{i}_dtemperature({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) {const_option} {{ return dthird_body_dtemperature({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature, {mixture_concentration});}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B, mixture_concentration = mixture_concentration)

def dthird_body_text_dlog_temperature(i, A, B, E, efficiencies, species_names, requires_mixture_concentration, configuration):
    mixture_concentration = "mixture_concentration"
    if requires_mixture_concentration:
        mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
    return_text = "{device_option}\n{scalar_function} dcall_forward_reaction_{i}_dlog_temperature({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) {const_option} {{ return dthird_body_dlog_temperature({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature, {mixture_concentration});}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B, mixture_concentration = mixture_concentration)

def dthird_body_text_dmixture_concentration(i, A, B, E, efficiencies, species_names, requires_mixture_concentration, dthird_body_multiplier_dspecies_index,configuration):
    mixture_concentration = "mixture_concentration"
    return_text = ''
    if requires_mixture_concentration:
        mixture_concentration = get_mixture_concentration(efficiencies, species_names, configuration)
        dmixture_concentration_dspecies = get_mixture_concentration_derivatives(efficiencies, species_names, configuration)
        return_text = "{device_option}\n{species_function} dcall_forward_reaction_{i}_dspecies({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) {const_option} {{ "+f"Species dmixture_concentration_dspecies = {{{{{dmixture_concentration_dspecies}}}}};" +"return scale_gen(dthird_body_dmixture_concentration({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature, {mixture_concentration}),dmixture_concentration_dspecies);}}"
    else:
        zero_species = ['{scalar_cast}(0)'.format(**vars(configuration))]*len(species_names)
        zero_species[dthird_body_multiplier_dspecies_index] = 'dthird_body_dmixture_concentration({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, log_temperature, {mixture_concentration})'
        single_species_return_text = '{{{joined}}}'.format(**vars(configuration), joined = ','.join(zero_species))
        return_text = "{device_option}\n{species_function} dcall_forward_reaction_{i}_dspecies({species_parameter} species, {scalar_parameter} temperature, {scalar_parameter} log_temperature, {scalar_parameter} mixture_concentration) {const_option} {{ return" + f" {{{single_species_return_text}}}" + ";}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B, mixture_concentration = mixture_concentration)


def create_reaction_functions_and_calls_third_body(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False, temperature_jacobian = False):
    if verbose:
        print(f"  Arrhenius Parameters (3-body reaction): A = {reaction.rate.pre_exponential_factor}, "
            f"b = {reaction.rate.temperature_exponent}, "
            f"Ea = {reaction.rate.activation_energy}")
        print(f"  Collision Partner Efficiencies: {reaction.efficiencies}")
    
    if reaction.third_body_name !='M':
        requires_mixture_concentration[reaction_index] = False
        third_body_index = species_names.index(reaction.third_body_name)
        dthird_body_multiplier_dspecies_index = third_body_index
        third_body_multiplier = "species[{index}]".format(index = third_body_index)
    else:
        requires_mixture_concentration[reaction_index] = True
        third_body_multiplier = 'mixture_concentration'
        dthird_body_multiplier_dspecies_index = None 
    reaction_rates[reaction_index] = third_body_text(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, reaction.efficiencies, species_names, requires_mixture_concentration[reaction_index], configuration)
    if temperature_jacobian:
        reaction_rates_derivatives.append(dthird_body_text_dtemperature(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, reaction.efficiencies, species_names, requires_mixture_concentration[reaction_index], configuration))
        reaction_rates_derivatives.append(dthird_body_text_dlog_temperature(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, reaction.efficiencies, species_names, requires_mixture_concentration[reaction_index], configuration))
    else:
        reaction_rates_derivatives.append(f'//dcall_forward_reaction_{reaction_index} temperature unused not needed')
    reaction_rates_derivatives.append(dthird_body_text_dmixture_concentration(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, reaction.efficiencies, species_names, requires_mixture_concentration[reaction_index], dthird_body_multiplier_dspecies_index, configuration))
    reaction_calls[reaction_index] = " call_forward_reaction_{reaction_index}(species, temperature, log_temperature, {third_body_multiplier});\n".format(**vars(configuration),reaction_index = reaction_index, third_body_multiplier = third_body_multiplier)

