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
        formatted_text = "{scalar_cast}({stoichiometric_production}) * rate_of_progress_{reaction_index}".format(**vars(configuration), stoichiometric_production = stoichiometric_production[index], reaction_index = reaction_index)
        if species_production_texts[index] == '':
            species_production_texts[index] = formatted_text
        else:
            species_production_texts[index] = ' + '.join([species_production_texts[index], formatted_text])

def get_reaction_function(reaction_rates, reaction_calls, reaction, configuration, reaction_index):
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

def create_rates_of_progress(progress_rates, reaction_index, forward_rate, backward_rate, configuration = None):
    if configuration == None:
        print("Warning this may cause compilation mismatch in decorators")
        configuration = get_configuration("configuration.yaml")
    formatted_text = (
    "{scalar} rate_of_progress_{reaction_index} = {forward_rate} * forward_reaction_{reaction_index} "
    ";//- {backward_rate} * forward_reaction_{reaction_index}/equilibrium_constant_{reaction_index};"
    .format(reaction_index=reaction_index, 
            forward_rate=forward_rate, 
            backward_rate=backward_rate, 
            **vars(configuration)))

    progress_rates[reaction_index] = formatted_text

def process_cantera_file(gas, configuration):
    species_names  = gas.species_names
    species_production_texts = [''] * gas.n_species
    species_production_jacobian = [[''] * (gas.n_species + 1)] * (gas.n_species + 1)
    reaction_rates = [''] * gas.n_reactions
    reaction_calls = [''] * gas.n_reactions
    progress_rates = [''] * gas.n_reactions
    [thermo_names, thermo_fits, thermo_types] = polyfit_thermodynamics(gas, configuration, order = int("{n_thermo_order}".format(**vars(configuration))))

    # Loop through all reactions
    for reaction_index in range(gas.n_reactions):
        reaction = gas.reaction(reaction_index)
        stoichiometric_production = np.zeros(len(species_names))
        stoichiometric_forward = np.zeros(len(species_names))
        stoichiometric_backward = np.zeros(len(species_names))
        indexes_of_species_in_reaction = []

        [forward_rate, backward_rate] = get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration = configuration)

        stoichiometric_production = stoichiometric_backward - stoichiometric_forward 

        
        accrue_species_production(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts, reaction_index)
        
        create_rates_of_progress(progress_rates, reaction_index, forward_rate, backward_rate, configuration)
        
        get_reaction_function(reaction_rates, reaction_calls, reaction, configuration, reaction_index)
    
    headers = []
    with open('types_inl.h','w') as file:
        write_type_defs(file, gas.n_species, configuration = configuration)
        headers.append('types_inl.h')

    with open('thermotransport_fits.h','w') as file:
        for name, thermo_fit, thermo_type in zip(thermo_names, thermo_fits, thermo_types):
            if thermo_type == "energy":
                write_energy_thermo_transport_fit(file, name, thermo_fit, configuration = configuration)
            else:
                write_thermo_transport_fit(file, name, thermo_fit,  configuration = configuration)
        headers.append('thermotransport_fits.h')
    
    with open('reactions.h','w') as file:
        write_reaction_rates(file, reaction_rates)
        headers.append('reactions.h')
    
    with open('source.h','w') as file:
        write_start_of_source_function(file, configuration=configuration)
        write_reaction_calculations(file, reaction_calls)
        write_progress_rates(file, progress_rates)
        write_species_production(file, species_production_texts, configuration = configuration)
        headers.append('source.h')
        write_end_of_function(file)
    
    required_headers = create_headers(configuration)
    return required_headers + headers


