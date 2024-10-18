import cantera as ct
import numpy as np
from .arrhenius import *
from .third_body import *
from .arithmetic import *
from .headers import *
from .configuration import *
from .thermo_chemistry import *
from .write import *
from .process_reactions import *

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

        [forward_rate, backward_rate] = get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration)
        stoichiometric_production = stoichiometric_backward - stoichiometric_forward 

        create_equilibrium_constants(stoichiometric_production, reaction_index, indexes_of_species_in_reaction, equilibrium_constants, configuration)
        accrue_species_production(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts, reaction_index, configuration)
        create_reaction_functions_and_calls(reaction_rates, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names)
        create_rates_of_progress(progress_rates, reaction_index, forward_rate, backward_rate, is_reversible, configuration)
    #sys.exit("Exiting the program")
    headers = []
    with open('types_inl.h','w') as file:
        write_type_defs(file, gas, configuration)
        headers.append('types_inl.h')

    with open('thermotransport_fits.h','w') as file:
        for name, thermo_fit, thermo_type in zip(thermo_names, thermo_fits, thermo_types):
            if thermo_type == "energy":
                write_energy_thermo_transport_fit(file, name, thermo_fit, configuration)
            elif thermo_type == "entropy":
                write_entropy_thermo_transport_fit(file, name, thermo_fit, configuration)
            elif thermo_type == "gibbs":
                write_gibbs_thermo_transport_fit(file, name, thermo_fit, configuration)
            else:
                write_thermo_transport_fit(file, name, thermo_fit,  configuration)
        headers.append('thermotransport_fits.h')
    
    with open('reactions.h','w') as file:
        write_reaction_rates(file, reaction_rates)
        headers.append('reactions.h')
    
    with open('source.h','w') as file:
        write_equilibrium_constants(file, equilibrium_constants, configuration)
        write_start_of_source_function(file, configuration=configuration)
        write_reaction_calculations(file, reaction_calls)
        write_progress_rates(file, progress_rates, is_reversible, equilibrium_constants, configuration)
        write_species_production(file, species_production_texts, configuration)
        headers.append('source.h')
        write_end_of_function(file)
    
    required_headers = create_headers(configuration)
    return required_headers + headers


