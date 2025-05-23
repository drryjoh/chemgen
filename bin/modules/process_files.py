import cantera as ct
import numpy as np
import importlib.util
from .arithmetic import *
from .headers import *
from .configuration import *
from .thermo_chemistry import *
from .write import *
from .process_reactions import *

def load_custom_sourcewriter(filepath):
    # Load a module from the given file path
    spec = importlib.util.spec_from_file_location("source_custom", filepath)
    custom_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(custom_module)
    # Return the SourceWriter class from the custom module
    return custom_module.SourceWriter

def process_cantera_file(gas, configuration, destination_folder, args, chemistry_solver, verbose = False, fit_gibbs_reaction = True, temperature_jacobian = False, remove_reactions = False):
    species_names  = gas.species_names
    species_production_texts = [''] * gas.n_species
    species_production_function_texts = [''] * gas.n_species
    species_production_on_fly_function_texts = [''] * gas.n_reactions
    species_production_jacobian_texts = [''] * gas.n_species
    species_production_jacobian_species_texts = [[] for _ in range(gas.n_species)]
    species_production_jacobian_species_indexes = [[] for _ in range(gas.n_species)]
    species_production_jacobian_temperature_texts = [''] * gas.n_species
    reaction_rates = [''] * gas.n_reactions
    reactions_depend_on = [[]] * gas.n_reactions
    reaction_rates_derivatives = []
    reaction_calls = [''] * gas.n_reactions
    progress_rates = [''] * gas.n_reactions
    progress_rates_derivatives = [''] * gas.n_reactions
    equilibrium_constants = [''] * gas.n_reactions
    dequilibrium_constants_dtemperature = [''] * gas.n_reactions
    is_reversible  = [False] * gas.n_reactions
    requires_mixture_concentration = [False] * gas.n_reactions  
    molecular_weights_string = ','.join(["{scalar_cast}({mw})".format(**vars(configuration), mw=mw) for mw in gas.molecular_weights])
    molecular_weights = f"{{{molecular_weights_string}}}"
    inv_molecular_weights_string = ','.join(["{scalar_cast}({mw})".format(**vars(configuration), mw=1.0/mw) for mw in gas.molecular_weights])
    inv_molecular_weights = f"{{{inv_molecular_weights_string}}}"
    [thermo_names, thermo_fits, thermo_types] = polyfit_thermodynamics(gas, configuration, order = int("{n_thermo_order}".format(**vars(configuration))))
    # Loop through all reactions
    for reaction_index in range(gas.n_reactions):
        reaction = gas.reaction(reaction_index)
        stoichiometric_production = np.zeros(len(species_names))
        stoichiometric_forward = np.zeros(len(species_names))
        stoichiometric_backward = np.zeros(len(species_names))
        indexes_of_species_in_reaction = []

        [forward_rate, backward_rate] = get_stoichmetric_balance_arithmetic(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration)
        [forward_rate_derivatives, backward_rate_derivatives] = get_stoichmetric_balance_arithmetic_derivatives(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration)
        stoichiometric_production = stoichiometric_backward - stoichiometric_forward 


        create_equilibrium_constants(stoichiometric_production, reaction_index, indexes_of_species_in_reaction, equilibrium_constants, dequilibrium_constants_dtemperature, configuration, fit_gibbs_reaction, temperature_jacobian = temperature_jacobian)

        accrue_species_production(indexes_of_species_in_reaction, stoichiometric_production, species_production_texts, species_production_function_texts, species_production_on_fly_function_texts, reaction_index, configuration)
        create_reaction_functions_and_calls(reaction_rates, reaction_rates_derivatives, reactions_depend_on, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose, temperature_jacobian = temperature_jacobian)
        create_rates_of_progress(progress_rates, species_production_function_texts, reaction_index, forward_rate, backward_rate, is_reversible, configuration) 
        create_rates_of_progress_derivatives(progress_rates_derivatives, reactions_depend_on, species_production_function_texts, reaction_index, forward_rate, backward_rate, forward_rate_derivatives, backward_rate_derivatives, is_reversible, indexes_of_species_in_reaction, stoichiometric_production, reaction, configuration, temperature_jacobian = temperature_jacobian)
        
        accrue_species_production_jacobian(indexes_of_species_in_reaction, stoichiometric_production, species_production_jacobian_species_texts,species_production_jacobian_species_indexes, species_production_jacobian_temperature_texts, reactions_depend_on, reaction_index, configuration, temperature_jacobian = temperature_jacobian)
    add_to_loops(species_production_jacobian_texts, species_production_jacobian_species_texts, species_production_jacobian_species_indexes, species_production_jacobian_temperature_texts, configuration, temperature_jacobian = False)
    headers = []
    with open(destination_folder/'types_inl.h','w') as file:
        write_type_defs(file, gas, configuration)
        headers.append('types_inl.h')
    with open(destination_folder/'generated_constants.h', 'w') as file:
        write_molecular_weights(file, molecular_weights, inv_molecular_weights,  configuration)
        write_species_names(file, species_names, configuration)
    with open(destination_folder/'thermotransport_fits.h','w') as file:
        for name, thermo_fit, thermo_type in zip(thermo_names, thermo_fits, thermo_types):
            if thermo_type == "energy":
                write_energy_thermo_transport_fit(file, name, thermo_fit, configuration)
            elif thermo_type == "entropy":
                write_entropy_thermo_transport_fit(file, name, thermo_fit, configuration)
            elif thermo_type == "gibbs":
                write_gibbs_thermo_transport_fit(file, name, thermo_fit, configuration)
            elif thermo_type == "gibbs_reaction":
                write_gibbs_reaction_transport_fit(file, name, thermo_fit, configuration)
            else:
                write_thermo_transport_fit(file, name, thermo_fit,  configuration)
    
    with open(destination_folder/'reactions.h','w') as file:
        write_reaction_rates(file, reaction_rates)
        write_reaction_rates_derivatives(file, reaction_rates_derivatives)
        headers.append('reactions.h')
    
    with open(destination_folder/'source.h','w') as file:
        if args.custom_source:
            try:
                # Load the custom SourceWriter
                CustomSourceWriter = load_custom_sourcewriter(args.custom_source)
                # Instantiate and use the custom SourceWriter
                custom_writer = CustomSourceWriter()
                custom_writer.write_source(file, equilibrium_constants, reaction_calls, 
                                        progress_rates, is_reversible, species_production_on_fly_function_texts,
                                        species_production_function_texts, 
                                        headers, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
                '''
                custom_writer.write_source_jacobian(file, equilibrium_constants, dequilibrium_constants_dtemperature, reactions_depend_on, reaction_calls, 
                                        progress_rates, progress_rates_derivatives, is_reversible, species_production_on_fly_function_texts,
                                        species_production_function_texts, 
                                        headers, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
                ## No Custom Writer For Jacobian
                '''    
            except (FileNotFoundError, AttributeError) as e:
                print(f"Error loading custom source writer: {e}")
                sys.exit(1)
        else:
            from .write_source import SourceWriter as source
            from .write_source_jacobian_species import SourceJacobianWriter as source_jacobian_species
            from .write_source_jacobian_remove_R import SourceJacobianWriter as source_jacobian_species_R
            from .write_source_jacobian_species_temperature import SourceJacobianWriter as source_jacobian_species_temperature
            source().write_source(file, equilibrium_constants, reaction_calls, progress_rates, is_reversible, species_production_on_fly_function_texts, species_production_texts, headers, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
            if temperature_jacobian:
                source_jacobian_species_temperature().write_source_jacobian(file, equilibrium_constants, dequilibrium_constants_dtemperature, reactions_depend_on,
                                                                            reaction_calls, progress_rates, progress_rates_derivatives, is_reversible, species_production_on_fly_function_texts, 
                                                                            species_production_texts, species_production_jacobian_texts, headers, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
            else:
                source_jacobian_species().write_source_jacobian(file, equilibrium_constants, dequilibrium_constants_dtemperature, reactions_depend_on,
                                                            reaction_calls, progress_rates, progress_rates_derivatives, is_reversible, species_production_on_fly_function_texts, 
                                                            species_production_texts, species_production_jacobian_texts, headers, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
                if remove_reactions:
                    source_jacobian_species_R().write_source_jacobian(file, equilibrium_constants, dequilibrium_constants_dtemperature, reactions_depend_on,
                                                                      reaction_calls, progress_rates, progress_rates_derivatives, is_reversible, species_production_on_fly_function_texts, 
                                                                      species_production_texts, species_production_jacobian_texts, headers, configuration, fit_gibbs_reaction =  fit_gibbs_reaction)
    
    required_headers = create_headers(configuration, chemistry_solver, destination_folder)
    
    return required_headers + headers


