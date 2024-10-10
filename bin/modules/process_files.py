import cantera as ct
import numpy as np
from .arrhenius import *
from .arithmetic import *
from .headers import *

def get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names):
    forward_rate_array = []
    for species, coeff in reaction.reactants.items():
        stoichiometric_forward[species_names.index(species)] = coeff
        forward_rate_array.append(raise_to_power(f"C{species_names.index(species)}", coeff))
        indexes_of_species_in_reaction.append(species_names.index(species))
    forward_rate = ' * '.join(forward_rate_array)

    backward_rate_array = []
    for species, coeff in reaction.products.items():
        stoichiometric_backward[species_names.index(species)] = coeff
        backward_rate_array.append(raise_to_power(f"C{species_names.index(species)}", coeff))
        indexes_of_species_in_reaction.append(species_names.index(species))
    backward_rate = ' * '.join(backward_rate_array)

    return (forward_rate, backward_rate)

def accrue_species_produciton(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts)
    for index in indexes_of_species_in_reaction: 
        if species_production_texts[index] == '':
            species_production_texts[index] = f"{stoichiometric_production[index]} * rate_of_progress_{i}"
        else:
            species_production_texts[index] = ' + '.join([species_production_texts[index], f"{stoichiometric_production[index]} * rate_of_progress_{i}"])

def get_reaction_function(reaction_rates, reaction, configuration):
        if isinstance(reaction, ct.Reaction):
            print(f"  Arrhenius Parameters: A = {reaction.rate.pre_exponential_factor}, "
                f"b = {reaction.rate.temperature_exponent}, "
                f"Ea = {reaction.rate.activation_energy}")
            reaction_rates[i] = arrhenius_text(i, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, configuration)
        
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

def process_cantera_file(gas, configuration = None):
    species_names  = gas.species_names
    species_production_texts = [''] * gas.n_species
    species_production_jacobian = [[''] * (gas.n_species + 1)] * (gas.n_species + 1)
    reaction_rates = [''] * gas.n_reactions
    progress_rates = [''] * gas.n_reactions

    # Loop through all reactions
    for i in range(gas.n_reactions):
        reaction = gas.reaction(i)
        stoichiometric_production = np.zeros(len(species_names))
        stoichiometric_forward = np.zeros(len(species_names))
        stoichiometric_backward = np.zeros(len(species_names))
        indexes_of_species_in_reaction = []

        # Print the reaction equation
        # print(f"Reaction {i + 1}: {reaction.equation}")

        [forward_rate, backward_rate] = get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names)

        stoichiometric_production = stoichiometric_backward - stoichiometric_forward 

        progress_rates[i] = f"rate_of_progress_{i} = {forward_rate} * r{i} - {backward_rate} * r{i}/equilibrium_constant_{i}"

        accrue_species_produciton(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts)

        get_reaction_function(reaction_rates, reaction, configuration)
        
        print()  # Add a blank line for better readability
    #write_reaction_rates(reaction_rates)
    #print(species_production_texts)
    #write_progress_rates(progress_rates)
    print(species_production_texts)
    clear_headers('./')
    create_headers(configuration = configuration)

