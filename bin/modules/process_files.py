import cantera as ct
import numpy as np
from .arrhenius import *
from .arithmetic import *

def process_cantera_file(gas, config=None):
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
        
        forward_rate = ''
        backward_rate = ''
        production_rate  = ''

        # Print the reaction equation
        print(f"Reaction {i + 1}: {reaction.equation}")

        print("reactants")
        forward_rate_array = []
        indexes_of_species_in_reaction = []
        for species, coeff in reaction.reactants.items():
            print(f"    index: {species_names.index(species)} {species}: {coeff}")
            stoichiometric_forward[species_names.index(species)] = coeff
            forward_rate_array.append(raise_to_power(f"C{species_names.index(species)}", coeff))
            indexes_of_species_in_reaction.append(species_names.index(species))

        forward_rate = ' * '.join(forward_rate_array)
        print(forward_rate)
        # Print stoichiometric coefficients for products
        print("  Products:")
        backward_rate_array = []
        for species, coeff in reaction.products.items():
            print(f"    {species}: {coeff}")
            stoichiometric_backward[species_names.index(species)] = coeff
            backward_rate_array.append(raise_to_power(f"C{species_names.index(species)}", coeff))
            indexes_of_species_in_reaction.append(species_names.index(species))
        backward_rate = ' * '.join(backward_rate_array)
        print(backward_rate)


        stoichiometric_production = stoichiometric_backward - stoichiometric_forward 
        progress_rates[i] = f"rate_of_progress_{i} = {forward_rate} * r{i} - {backward_rate} * r{i}/equilibrium_constant_{i}"
        for index in indexes_of_species_in_reaction: 
            if species_production_texts[index] == '':
                species_production_texts[index] = f"{stoichiometric_production[index]} * rate_of_progress_{i}"
            else:
                species_production_texts[index] = ' + '.join([species_production_texts[index], f"{stoichiometric_production[index]} * rate_of_progress_{i}"])

        print(f"stoichiometric_production: {stoichiometric_production}")
        print(f"stoichiometric_forward: {stoichiometric_forward}")
        print(f"stoichiometric_backward: {stoichiometric_backward}")

        # Get the type of the reaction
        reaction_type = type(reaction).__name__
        print(f"Reaction Type: {reaction_type}")
        
        print(reaction_type)
        # Print different parameters depending on the reaction type
        if isinstance(reaction, ct.Reaction):
            print(f"  Arrhenius Parameters: A = {reaction.rate.pre_exponential_factor}, "
                f"b = {reaction.rate.temperature_exponent}, "
                f"Ea = {reaction.rate.activation_energy}")
            reaction_rates[i] = arrhenius_text(i, reaction.rate.pre_exponential_factor, reaction.rate.temperature_exponent, reaction.rate.activation_energy, config)
        
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
        
        print()  # Add a blank line for better readability
    #write_reaction_rates(reaction_rates)
    #print(species_production_texts)
    #write_progress_rates(progress_rates)
    print(species_production_texts)
