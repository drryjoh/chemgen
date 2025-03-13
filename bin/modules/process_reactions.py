from .reactions.reactions import *
from .arithmetic import *
from .headers import *
from .configuration import *
from .thermo_chemistry import *
from .write import *
import sys

def get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration):
    forward_rate_array = []
    for species, coeff in reaction.reactants.items():
        stoichiometric_forward[species_names.index(species)] = coeff
        species_index = species_names.index(species)
        species_element_i  = configuration.species_element.format(i = species_index)
        forward_rate_array.append(raise_to_power(species_element_i, coeff))
        indexes_of_species_in_reaction.append(species_names.index(species))
    forward_rate = ' * '.join(forward_rate_array)

    backward_rate_array = []
    for species, coeff in reaction.products.items():
        stoichiometric_backward[species_names.index(species)] = coeff
        species_index = species_names.index(species)
        species_element_i  = configuration.species_element.format(i = species_index)
        backward_rate_array.append(raise_to_power(species_element_i, coeff))
        indexes_of_species_in_reaction.append(species_names.index(species))
    backward_rate = ' * '.join(backward_rate_array)

    return (forward_rate, backward_rate)
#species_production_on_fly_function_texts,
def accrue_species_production(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts, species_production_function_texts,  species_production_on_fly_function_texts, reaction_index, configuration):
    for index in indexes_of_species_in_reaction: 
        formatted_text = "{scalar_cast}({stoichiometric_production}) * rate_of_progress_{reaction_index}".format(**vars(configuration), 
        stoichiometric_production = stoichiometric_production[index], 
        reaction_index = reaction_index)
        if species_production_texts[index] == '':
            species_production_texts[index] = formatted_text
        else:
            species_production_texts[index] = ' + '.join([species_production_texts[index], formatted_text])
    formatted_text = ''

    for index in indexes_of_species_in_reaction: 
        formatted_text = "{scalar_cast}({stoichiometric_production}) * progress_rates[{reaction_index}]".format(**vars(configuration), 
        stoichiometric_production = stoichiometric_production[index], 
        reaction_index = reaction_index)
        if species_production_function_texts[index] == '':
            species_production_function_texts[index] = formatted_text
        else:
            species_production_function_texts[index] = ' + '.join([species_production_function_texts[index], formatted_text])
    
    on_the_fly_production = []
    for index in indexes_of_species_in_reaction: 
        on_the_fly_production.append( "species_source[{species_index}] = {scalar_cast}({stoichiometric_production}) * progress_rate; //Reaction {reaction_index}".format(**vars(configuration), 
        stoichiometric_production = stoichiometric_production[index], 
        reaction_index = reaction_index,
        species_index = index))
    species_production_on_fly_function_texts[reaction_index] = '\n'.join(on_the_fly_production)

def create_reaction_functions_and_calls(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False):
    is_reversible[reaction_index] = reaction.reversible
    if reaction.reaction_type == "Arrhenius":
        create_reaction_functions_and_calls_arrhenius(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose)
    elif reaction.reaction_type == "three-body-Arrhenius":
        create_reaction_functions_and_calls_third_body(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose)
    elif "falloff" in reaction.reaction_type:
        create_reaction_functions_and_calls_falloff(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose)
    elif reaction.reaction_type == "pressure-dependent-Arrhenius":
        create_reaction_functions_and_calls_pressure_dependent_arrhenius(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose)
    elif reaction.reaction_type == "Chebyshev":
        print(f"  Chebyshev Reaction Coefficients:")
        print(f"    Tmin = {reaction.Tmin}, Tmax = {reaction.Tmax}")
        print(f"    Pmin = {reaction.Pmin}, Pmax = {reaction.Pmax}")
        print(f"    Coefficients: {reaction.coeffs}")

    else:
        print(f"  Unknown reaction type: {reaction.reaction_type }")

def create_rates_of_progress(progress_rates, progress_rates_functions, reaction_index, forward_rate, backward_rate, is_reversible, configuration):
    if is_reversible[reaction_index]:
        formatted_text = (
        "{scalar} rate_of_progress_{reaction_index} = {forward_rate} * forward_reaction_{reaction_index} "
        "- {backward_rate} * forward_reaction_{reaction_index}/equilibrium_constant_{reaction_index};"
        .format(reaction_index=reaction_index, 
                forward_rate=forward_rate, 
                backward_rate=backward_rate, 
                **vars(configuration)))
    else:
        formatted_text = (
        "{scalar} rate_of_progress_{reaction_index} = {forward_rate} * forward_reaction_{reaction_index};"
        .format(reaction_index=reaction_index, 
                forward_rate=forward_rate,
                **vars(configuration)))

    progress_rates[reaction_index] = formatted_text
def create_equilibrium_constants(stoichiometric_production, reaction_index, indexes_of_species_in_reaction, equilibrium_constants, dequilibrium_constants_dtemperature, configuration, fit_gibbs_reaction = True):
    scalar_cast = "{scalar_cast}".format(**vars(configuration))
    equilibrium_constant_elements = []
    sum_stoichiometric_production = np.sum(stoichiometric_production)
    
    power_term = ''
    dpower_term_dtemperature = ''
    if sum_stoichiometric_production.is_integer():
        power_integer = int(sum_stoichiometric_production)
        if power_integer < 0:
            power_term = raise_to_power('inv_pressure_atmosphere() * universal_gas_constant() * temperature', np.abs(power_integer))
            dpower_term_dtemperature = draise_to_power_chain('inv_pressure_atmosphere() * universal_gas_constant() * temperature', 'inv_pressure_atmosphere() * universal_gas_constant()', np.abs(power_integer))
        elif power_integer > 0:
            power_term = raise_to_power("pressure_atmosphere() * inv_universal_gas_constant_temperature",power_integer)
            dpower_term_dtemperature = draise_to_power_chain("pressure_atmosphere() * inv_universal_gas_constant_temperature", "pressure_atmosphere() * dinv_universal_gas_constant_temperature_dtemperature", np.abs(power_integer))
        elif power_integer == 0:
            power_term = f"{scalar_cast}(1.0)"
            dpower_term_dtemperature = f"{scalar_cast}(0.0)"
        else:
            power_term = f'pow_gen(pressure_atmosphere() * inv_universal_gas_constant_temperature, {power_integer})'
    else:
        power_term = f'pow_gen(pressure_atmosphere() * inv_universal_gas_constant_temperature, {scalar_cast}({sum_stoichiometric_production}))'
        dpower_term_dtemperature = f'multiply(dpow_gen_da(pressure_atmosphere() * inv_universal_gas_constant_temperature, {scalar_cast}({sum_stoichiometric_production})), pressure_atmosphere() * dinv_universal_gas_constant_temperature_dtemperature)'
    
    for index in indexes_of_species_in_reaction:
        if stoichiometric_production[index] != 0:
            equilibrium_constant_elements.append(f"{scalar_cast}({stoichiometric_production[index]}) * gibbs_free_energies[{index}]")

    if fit_gibbs_reaction:
        equilibrium_constants[reaction_index] = "multiply(exp_gen(-gibbs_reactions[{reaction_index}]), {power_term})".format(power_term = power_term, reaction_index = reaction_index)
        dequilibrium_constants_dtemperature[reaction_index] = "multiply_chain(exp_gen(-gibbs_reactions[{reaction_index}]), exp_chain(-gibbs_reactions[{reaction_index}], -dgibbs_reactions_dlog_temperature[{reaction_index}]) * dlog_temperature_dtemperature, {power_term}, {dpower_term_dtemperature})".format(power_term = power_term, reaction_index = reaction_index)
    else:
        equilibrium_constants[reaction_index] = "exp_gen(-({gibbs_sum}) * inv_universal_gas_constant_temperature) * {power_term}".format(gibbs_sum = '+'.join(equilibrium_constant_elements).replace("+-","-"), power_term=power_term)
