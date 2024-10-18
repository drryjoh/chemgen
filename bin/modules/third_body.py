def third_body_text(i, A, B, E, efficiencies, species_names, configuration):
    mixture_concentration_array = []

    if not efficiencies:
        efficiencies = {specie: 1.0 for specie in species_names}

    for species_index, specie in enumerate(species_names):
        if specie in efficiencies:
            if efficiencies[specie] != 0:
                mixture_concentration_array.append(f"{configuration.scalar_cast}({efficiencies[specie]})*{configuration.species_element.format(i = species_index)}")
        mixture_concentration_array.append(f"{configuration.scalar_cast}(1)*{configuration.species_element.format(i = species_index)}")

    return_text = "{device_option}\n{scalar_function} call_forward_reaction_{i}({species_parameter} species, {scalar_parameter} temperature) {const_option} {{ return third_body({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature, {mixture_concentration});}}"
    return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B, mixture_concentration = ' + '.join(mixture_concentration_array))

def dthird_body_dtemperature_text(i, A, B, E, configuration):
    return f"auto dreaction_dtemperature_{i}(const double& temperature) const {{ return dthird_body_dtemperature({A}, {B}, {E}, temperature)}}"

def create_reaction_functions_and_calls_third_body(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names):
    print(f"  Arrhenius Parameters (3-body reaction): A = {reaction.rate.pre_exponential_factor}, "
        f"b = {reaction.rate.temperature_exponent}, "
        f"Ea = {reaction.rate.activation_energy}")
    print(f"  Collision Partner Efficiencies: {reaction.efficiencies}")
    requires_mixture_concentration[reaction_index] = True
    reaction_rates[reaction_index] = third_body_text(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, reaction.efficiencies, species_names, configuration)
    reaction_calls[reaction_index] = "{scalar} forward_reaction_{reaction_index} = call_forward_reaction_{reaction_index}(species, temperature);\n".format(**vars(configuration),reaction_index = reaction_index)

