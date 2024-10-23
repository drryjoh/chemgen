from .falloff_troe import *
from .falloff_lindemann import *
import sys

def create_reaction_functions_and_calls_falloff(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names):
    falloff_type = reaction.reaction_type
    if "Troe" in falloff_type:
        create_reaction_functions_and_calls_troe(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names)
    elif "Lindemann" in falloff_type:
        if reaction.rate.falloff_coeffs>0:
            print(f"Expected coefficients for {falloff_type} is incorrect")
            sys.exit()
        else:
            print("good job")
            create_reaction_functions_and_calls_lindemann(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names)
    else:
        print(f"Support for {falloff_type} is missing")
        print(reaction.rate.falloff_coeffs)
        sys.exit()
    '''

    elif falloff_type == "SRI":
        print(f"  This is an SRI falloff reaction with parameters: {reaction.falloff.parameters}")
    
    elif falloff_type == "Lindemann":
        print(f"  This is a Lindemann falloff reaction (no additional parameters).")
    
    else:
        print(f"  Unknown falloff model: {falloff_type}")
    '''
