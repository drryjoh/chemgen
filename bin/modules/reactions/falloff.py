from .falloff_troe import *
from .falloff_lindemann import *
from .falloff_sri import *
import sys

def create_reaction_functions_and_calls_falloff(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False):
    falloff_type = reaction.reaction_type
    if "Troe" in falloff_type:
        create_reaction_functions_and_calls_troe(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose)
    elif "Lindemann" in falloff_type:
        if reaction.rate.falloff_coeffs>0:
            print(f"Expected coefficients for {falloff_type} is incorrect")
        else:
            create_reaction_functions_and_calls_lindemann(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose)
    elif "SRI" in falloff_type:
        create_reaction_functions_and_calls_sri(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose)
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
