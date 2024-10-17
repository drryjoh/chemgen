import cantera as ct
import numpy as np
from .arrhenius import *
from .arithmetic import *
from .headers import *
from .configuration import *
from .thermo_chemistry import *
from .write import *

def get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration = None):
    if configuration == None:
        configuration = get_configuration("configuration.yaml")
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

def accrue_species_production(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts, reaction_index, configuration = None):
    if configuration == None:
        print("Warning this may cause compilation mismatch in decorators")
        configuration = get_configuration("configuration.yaml")

    for index in indexes_of_species_in_reaction: 
        formatted_text = "{scalar_cast}({stoichiometric_production}) * rate_of_progress_{reaction_index}".format(**vars(configuration), 
        stoichiometric_production = stoichiometric_production[index], 
        reaction_index = reaction_index)
        if species_production_texts[index] == '':
            species_production_texts[index] = formatted_text
        else:
            species_production_texts[index] = ' + '.join([species_production_texts[index], formatted_text])

def get_reaction_function(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration):
        is_reversible[reaction_index] = reaction.reversible
        if isinstance(reaction, ct.Reaction):
            print(f"  Arrhenius Parameters: A = {reaction.rate.pre_exponential_factor}, "
                f"b = {reaction.rate.temperature_exponent}, "
                f"Ea = {reaction.rate.activation_energy}")
            reaction_rates[reaction_index] = arrhenius_text(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, configuration)
            reaction_calls[reaction_index] = "{scalar} forward_reaction_{reaction_index} = call_forward_reaction_{reaction_index}(temperature);\n".format(**vars(configuration),reaction_index = reaction_index)
        
        elif isinstance(reaction, ct.ThreeBodyReaction):
            print(f"  Arrhenius Parameters (3-body reaction): A = {reaction.rate.pre_exponential_factor}, "
                f"b = {reaction.rate.temperature_exponent}, "
                f"Ea = {reaction.rate.activation_energy}")
            print(f"  Collision Partner Efficiencies: {reaction.efficiencies}")
            requires_mixture_concentration[reaction_index] = True
            reaction_rates[reaction_index] = three_body_text(reaction_index, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, configuration)
        
        elif isinstance(reaction, ct.FalloffReaction):
            print(f"  Arrhenius Parameters (high pressure limit): A = {reaction.high_rate.pre_exponential_factor}, "
                f"b = {reaction.high_rate.temperature_exponent}, "
                f"Ea = {reaction.high_rate.activation_energy}")
            print(f"  Arrhenius Parameters (low pressure limit): A = {reaction.low_rate.pre_exponential_factor}, "
                f"b = {reaction.low_rate.temperature_exponent}, "
                f"Ea = {reaction.low_rate.activation_energy}")
            print(f"  Falloff Parameters: {reaction.falloff.parameters}")
        
        elif isinstance(reaction, ct.PlogReaction):
            print("  PLOG Reaction with Rate Expressions:")
            for P, rate in reaction.rates:
                print(f"    At {P} Pa: A = {rate.pre_exponential_factor}, "
                    f"b = {rate.temperature_exponent}, "
                    f"Ea = {rate.activation_energy}")
        
        elif isinstance(reaction, ct.ChebyshevReaction):
            print(f"  Chebyshev Reaction Coefficients:")
            print(f"    Tmin = {reaction.Tmin}, Tmax = {reaction.Tmax}")
            print(f"    Pmin = {reaction.Pmin}, Pmax = {reaction.Pmax}")
            print(f"    Coefficients: {reaction.coeffs}")

        else:
            print(f"  Unknown reaction type: {reaction_type}")

def create_rates_of_progress(progress_rates, reaction_index, forward_rate, backward_rate, is_reversible, configuration):
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
def create_equilibrium_constants(stoichiometric_production, reaction_index, indexes_of_species_in_reaction, equilibrium_constants, configuration):
    scalar_cast = "{scalar_cast}".format(**vars(configuration))
    equilibrium_constant_elements = []
    sum_stoichiometric_production = np.sum(stoichiometric_production)
    
    #pow(p_atm*inv(R*T),{power_integer})
    power_term = ''
    print(sum_stoichiometric_production.is_integer())
    if sum_stoichiometric_production.is_integer():
        power_integer = int(sum_stoichiometric_production)
        print(power_integer)
        if power_integer < 0:
            power_term = raise_to_power('inv_pressure_atmosphere() * universal_gas_constant() * temperature', np.abs(power_integer))
        elif power_integer > 0:
            power_term = raise_to_power("pressure_atmosphere() * inv_universal_gas_constant_temperature",power_integer)
        elif power_integer == 0:
            power_term = f"{scalar_cast}(1.0)"
        else:
            power_term = f'pow_gen(pressure_atmosphere() * inv_universal_gas_constant_temperature,{power_integer})'
    else:
        power_term = f'pow_gen(pressure_atmosphere() * inv_universal_gas_constant_temperature,{scalar_cast}({sum_stoichiometric_production}))'
    print(power_term)
    for index in indexes_of_species_in_reaction:
        equilibrium_constant_elements.append(f"{scalar_cast}({stoichiometric_production[index]}) * gibbs_free_energies[{index}]")
    equilibrium_constants[reaction_index] = "exp_gen(-({gibbs_sum}) * inv_universal_gas_constant_temperature) * {power_term}".format(gibbs_sum = '+'.join(equilibrium_constant_elements).replace("+-","-"),
    power_term=power_term)

def process_cantera_file(gas, configuration):
    species_names  = gas.species_names
    species_production_texts = [''] * gas.n_species
    species_production_jacobian = [[''] * (gas.n_species + 1)] * (gas.n_species + 1)
    reaction_rates = [''] * gas.n_reactions
    reaction_calls = [''] * gas.n_reactions
    progress_rates = [''] * gas.n_reactions
    equilibrium_constants = [''] * gas.n_reactions
    is_reversible  = [False] * gas.n_reactions
    requires_mixture_concentration = [False] * gas.n_reactions

    [thermo_names, thermo_fits, thermo_types] = polyfit_thermodynamics(gas, configuration, order = int("{n_thermo_order}".format(**vars(configuration))))

    # Loop through all reactions
    for reaction_index in range(gas.n_reactions):
        reaction = gas.reaction(reaction_index)
        stoichiometric_production = np.zeros(len(species_names))
        stoichiometric_forward = np.zeros(len(species_names))
        stoichiometric_backward = np.zeros(len(species_names))
        indexes_of_species_in_reaction = []

        [forward_rate, backward_rate] = get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration = configuration)
        print(f"{reaction_index} here")
        print(stoichiometric_backward)
        print(stoichiometric_forward)

        stoichiometric_production = stoichiometric_backward - stoichiometric_forward 

        create_equilibrium_constants(stoichiometric_production, reaction_index, indexes_of_species_in_reaction, equilibrium_constants, configuration)

        accrue_species_production(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts, reaction_index)
        get_reaction_function(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration)
        create_rates_of_progress(progress_rates, reaction_index, forward_rate, backward_rate, is_reversible, configuration)
    
    headers = []
    with open('types_inl.h','w') as file:
        write_type_defs(file, gas, configuration = configuration)
        headers.append('types_inl.h')

    with open('thermotransport_fits.h','w') as file:
        for name, thermo_fit, thermo_type in zip(thermo_names, thermo_fits, thermo_types):
            if thermo_type == "energy":
                write_energy_thermo_transport_fit(file, name, thermo_fit, configuration = configuration)
            elif thermo_type == "entropy":
                write_entropy_thermo_transport_fit(file, name, thermo_fit, configuration = configuration)
            elif thermo_type == "gibbs":
                write_gibbs_thermo_transport_fit(file, name, thermo_fit, configuration = configuration)
            else:
                write_thermo_transport_fit(file, name, thermo_fit,  configuration = configuration)
        headers.append('thermotransport_fits.h')
    
    with open('reactions.h','w') as file:
        write_reaction_rates(file, reaction_rates)
        headers.append('reactions.h')
    
    with open('source.h','w') as file:
        write_equilibrium_constants(file, equilibrium_constants, configuration)
        write_start_of_source_function(file, configuration=configuration)
        write_reaction_calculations(file, reaction_calls)
        write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration)
        write_species_production(file, species_production_texts, configuration = configuration)
        headers.append('source.h')
        write_end_of_function(file)
    
    required_headers = create_headers(configuration)
    return required_headers + headers


