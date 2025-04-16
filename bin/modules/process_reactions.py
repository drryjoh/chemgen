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

def get_stoichmetric_balance_arithmetic_derivatives(stoichiometric_forward, stoichiometric_backward, indexes_of_species_in_reaction, reaction, species_names, configuration):
    
    forward_rate_array = []
    dforward_rate_array = []
    forward_indices = []
    for species, coeff in reaction.reactants.items():
        species_index = species_names.index(species)
        forward_indices.append(species_index)
        species_element_i  = configuration.species_element.format(i = species_index)
        forward_rate_array.append(raise_to_power(species_element_i, coeff))

    backward_rate_array = []
    dbackward_rate_array = []
    backward_indices = []
    for species, coeff in reaction.products.items():
        species_index = species_names.index(species)
        backward_indices.append(species_index)
        species_element_i  = configuration.species_element.format(i = species_index)
        backward_rate_array.append(raise_to_power(species_element_i, coeff))
    
    for species, coeff in reaction.reactants.items():
        species_index = species_names.index(species)
        species_element_i  = configuration.species_element.format(i = species_index)
        dforward_rate_array.append(draise_to_power(species_element_i, coeff))

    for species, coeff in reaction.products.items():
        species_index = species_names.index(species)
        species_element_i  = configuration.species_element.format(i = species_index)
        dbackward_rate_array.append(draise_to_power(species_element_i, coeff))
    
    dforward_rates = ['0']*len(species_names)
    dbackward_rates = ['0']*len(species_names)

    for i, forward_index in enumerate(forward_indices):
        other_species = forward_rate_array.copy()
        other_species.pop(i)
        if len(other_species) > 1:
            dforward_rates[forward_index] = '{0} * {1}'.format(dforward_rate_array[i], '*'.join(other_species))
        elif len(other_species) == 1:
            dforward_rates[forward_index] = '{0} * {1}'.format(dforward_rate_array[i], other_species[0])
        else:
            dforward_rates[forward_index] = '{0}'.format(dforward_rate_array[i])

    for i, backward_index in enumerate(backward_indices):
        other_species = backward_rate_array.copy()
        other_species.pop(i)
        if len(other_species) > 1:
            dbackward_rates[backward_index] = '{0} * {1}'.format(dbackward_rate_array[i], '*'.join(other_species))
        elif len(other_species) == 1:
            dbackward_rates[backward_index] = '{0} * {1}'.format(dbackward_rate_array[i], other_species[0])
        else:
            dbackward_rates[backward_index] = '{0}'.format(dbackward_rate_array[i])
    for i, db_rates in enumerate(dbackward_rates):
        dbackward_rates[i].replace('1 * ','')
        dforward_rates[i].replace('1 * ','')

    return (dforward_rates, dbackward_rates)


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

def accrue_species_production_jacobian(indexes_of_species_in_reaction, stoichiometric_production, species_production_jacobian_species_texts, species_production_jacobian_species_indexes, species_production_jacobian_temperature_texts, reactions_depend_on, reaction_index, configuration, temperature_jacobian = False):
    depends_on_temperature = True; depends_on_species = True
    begin = '0'
    end = 'n_species'
    if temperature_jacobian:
        begin = '1'
        end = 'n_species + 1'
    for index in indexes_of_species_in_reaction:
        if depends_on_temperature:
            if temperature_jacobian:
                if species_production_jacobian_temperature_texts[index] == "":
                    species_production_jacobian_temperature_texts[index] = "{scalar_cast}({stoichiometric_production}) * drate_of_progress_{reaction_index}_dtemperature".format(**vars(configuration), stoichiometric_production = stoichiometric_production[index], reaction_index = reaction_index)  
                else:
                    species_production_jacobian_temperature_texts[index] += "+ {scalar_cast}({stoichiometric_production}) * drate_of_progress_{reaction_index}_dtemperature".format(**vars(configuration), stoichiometric_production = stoichiometric_production[index], reaction_index = reaction_index)  
        if depends_on_species:
                species_production_jacobian_species_texts[index].append( "{stoichiometric_production}".format(**vars(configuration), stoichiometric_production = stoichiometric_production[index], reaction_index = reaction_index)) 
                species_production_jacobian_species_indexes[index].append("{reaction_index}".format(reaction_index= reaction_index))


def add_to_loops(species_production_jacobian_texts, species_production_jacobian_species_texts, species_production_jacobian_species_indexes, species_production_jacobian_temperature_texts, configuration, temperature_jacobian = False):
    begin = '0'
    end = 'n_species'
    jacobian_temperature = ""
    if temperature_jacobian:
        begin = '1'
        end = 'n_species + 1'
    
    for i, jacobian_texts in enumerate(species_production_jacobian_species_texts):
        if temperature_jacobian:
            jacobian_temperature = "jacobian_net_production_rates[{species_index}+{begin}][0] = {jacobian_temperature_text};".format(species_index = i, begin =begin, jacobian_temperature_text = species_production_jacobian_temperature_texts[i])
        else:
            jacobian_temperature  = "        //no temperature jacobian"
        if jacobian_texts != []:
            n_coefficients = len(jacobian_texts)
            jacobian_species_texts = """
            const {index} n_rates_of_progres_species_jacobian_{species_index} = {n_coefficients};
            static constexpr {scalar_list}<{scalar}, n_rates_of_progres_species_jacobian_{species_index}> coefficients_{species_index} = {{{coeffs}}};
            static constexpr {scalar_list}<{index}, n_rates_of_progres_species_jacobian_{species_index}> idx_{species_index} = {{{indexes}}};
            """.format(**vars(configuration), n_coefficients = n_coefficients, coeffs = ','.join(jacobian_texts), indexes = ",".join(species_production_jacobian_species_indexes[i]), species_index = i)
            species_production_jacobian_texts[i] = """
        /*
        {jacobian_temperature}
        {jacobian_species_texts}
        for({index} i = {begin}; i < {end}; i++)
        {{
            {scalar} sum = 0;
            for({index} j = 0; j < n_rates_of_progres_species_jacobian_{species_index}; j++)
            {{
                sum += coefficients_{species_index}[j] * drate_of_progress_dspecies[idx_{species_index}[j]][i];
            }}
            jacobian_net_production_rates[{species_index}][i] = sum;
        }}
        */
            """.format(**vars(configuration), species_index = i, begin = begin, end = end, jacobian_texts = jacobian_texts, jacobian_temperature = jacobian_temperature, jacobian_species_texts = jacobian_species_texts)  
        else:
            species_production_jacobian_texts[i] ="""
        //no species jacobian
            """
def create_reaction_functions_and_calls(reaction_rates, reaction_rates_derivatives, reactions_depend_on, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = False, temperature_jacobian = False):
    is_reversible[reaction_index] = reaction.reversible
    
    if reaction.reaction_type == "Arrhenius":
        reactions_depend_on[reaction_index] = ['temperature','log_temperature']
        create_reaction_functions_and_calls_arrhenius(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose, temperature_jacobian = temperature_jacobian)
    elif reaction.reaction_type == "three-body-Arrhenius":
        reactions_depend_on[reaction_index] = ['temperature','species','log_temperature']
        create_reaction_functions_and_calls_third_body(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose, temperature_jacobian = temperature_jacobian)
    elif "falloff" in reaction.reaction_type:
        reactions_depend_on[reaction_index] = ['temperature','species','log_temperature']
        create_reaction_functions_and_calls_falloff(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose, temperature_jacobian = temperature_jacobian)
    elif reaction.reaction_type == "pressure-dependent-Arrhenius":
        reactions_depend_on[reaction_index] = ['temperature','pressure']
        create_reaction_functions_and_calls_pressure_dependent_arrhenius(reaction_rates, reaction_rates_derivatives, reaction_calls, reaction, configuration, reaction_index, is_reversible, requires_mixture_concentration, species_names, verbose = verbose, temperature_jacobian = temperature_jacobian)
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
        "{scalar} rate_of_progress_{reaction_index} = multiply({forward_rate}, forward_reaction_{reaction_index}) "
        "- multiply({backward_rate}, divide(forward_reaction_{reaction_index}, equilibrium_constant_{reaction_index}));"
        .format(reaction_index = reaction_index, 
                forward_rate = forward_rate, 
                backward_rate = backward_rate, 
                **vars(configuration)))
    else:
        formatted_text = (
        "{scalar} rate_of_progress_{reaction_index} = multiply({forward_rate}, forward_reaction_{reaction_index});"
        .format(reaction_index = reaction_index, 
                forward_rate = forward_rate,
                **vars(configuration)))

    progress_rates[reaction_index] = formatted_text
def add_to_jacobian(variable, index_with_respect_to, indexes_of_species_in_reaction, stoichiometric_production):
    running_text = []
    for i, species_index in enumerate(indexes_of_species_in_reaction):
        if stoichiometric_production[species_index]!=0:
            running_text.append(f"        jacobian_net_production_rates[{species_index}][{index_with_respect_to}] += {stoichiometric_production[species_index]}*{variable};\n")
    return ''.join(running_text)
def add_to_jacobian_all(variable, indexes_of_species_in_reaction, stoichiometric_production):
    running_text = []
    for i, species_index in enumerate(indexes_of_species_in_reaction):
        if stoichiometric_production[species_index]!=0:
            running_text.append(f"        jacobian_net_production_rates[{species_index}] = jacobian_net_production_rates[{species_index}] + scale_gen({stoichiometric_production[species_index]}, {variable});\n")
    return ''.join(running_text)
def create_rates_of_progress_derivatives(progress_rates_derivatives, reactions_depend_on, progress_rates_functions, 
                                         reaction_index, forward_rate, backward_rate, forward_rate_derivatives,
                                         backward_rate_derivatives, is_reversible, 
                                         indexes_of_species_in_reaction, stoichiometric_production, reaction,
                                         configuration, temperature_jacobian = False):
    formatted_text = ''
    if is_reversible[reaction_index]:
        if temperature_jacobian:
            if "temperature" in reactions_depend_on[reaction_index]:
                formatted_text += """
        {scalar} drate_of_progress_{reaction_index}_dtemperature = 
        multiply({forward_rate}, dforward_reaction_{reaction_index}_dtemperature)
        - 
        multiply({backward_rate}, 
                 divide_chain(forward_reaction_{reaction_index},
                              dforward_reaction_{reaction_index}_dtemperature, 
                              equilibrium_constant,
                              dequilibrium_constant_dtemperature));
""".format(reaction_index = reaction_index, 
                    forward_rate = forward_rate, 
                    backward_rate = backward_rate, 
                    **vars(configuration))
        else:
            formatted_text += """
        //drate_of_progress_temperature unused
            """
        if ("species" in reactions_depend_on[reaction_index]) or ("pressure" in reactions_depend_on[reaction_index]):
            dforward_rate_dspecies = ''
            dbackward_rate_dspecies = ''
            for species_index, dforward_rate in enumerate(forward_rate_derivatives):
                #jacobian_net_production_rates
                if dforward_rate == '1':
                    dforward_rate_dspecies += f"        drate_of_progress_dspecies = forward_reaction_{reaction_index};// {reaction_index} {species_index} \n"
                    dforward_rate_dspecies += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
                elif dforward_rate != '0':
                    dforward_rate_dspecies += f"        drate_of_progress_dspecies = multiply({dforward_rate}, forward_reaction_{reaction_index});//{reaction_index} {species_index}\n"
                    dforward_rate_dspecies += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)

            for species_index, dbackward_rate in enumerate(backward_rate_derivatives):
                if dbackward_rate == '1':
                    dforward_rate_dspecies += f"        drate_of_progress_dspecies = -divide(forward_reaction_{reaction_index}, equilibrium_constant); //{reaction_index} {species_index}\n"
                    dforward_rate_dspecies += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
                elif dbackward_rate != '0':
                    dforward_rate_dspecies += f"        drate_of_progress_dspecies = -multiply({dbackward_rate}, divide(forward_reaction_{reaction_index}, equilibrium_constant));// {reaction_index} {species_index}\n"
                    dforward_rate_dspecies += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
            
            formatted_text += """
{dforward_rate_dspecies}
{dbackward_rate_dspecies}
        drate_of_progress_dspecies_all_species = 
        scale_gen({forward_rate}, dforward_reaction_{reaction_index}_dspecies) -
        scale_gen(divide({backward_rate}, 
                         equilibrium_constant), 
                  dforward_reaction_{reaction_index}_dspecies);
{all_species}
                        """.format(reaction_index = reaction_index, 
                    forward_rate = forward_rate, 
                    backward_rate = backward_rate, 
                    dforward_rate_dspecies = dforward_rate_dspecies,
                    dbackward_rate_dspecies = dbackward_rate_dspecies,
                    all_species = add_to_jacobian_all("drate_of_progress_dspecies_all_species", indexes_of_species_in_reaction, stoichiometric_production),
                    **vars(configuration))
        else:
            formatted_text += """        //drate_of_progress_dspecies[{reaction_index}] = {{{scalar_cast}(0)}};
""".format(reaction_index = reaction_index, **vars(configuration))
            for species_index, dforward_rate in enumerate(forward_rate_derivatives):
                if dforward_rate == '1':
                    formatted_text += f"        drate_of_progress_dspecies = forward_reaction_{reaction_index}; // [{reaction_index}][{species_index}]\n"
                    formatted_text += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
                elif dforward_rate != '0':
                    formatted_text += f"        drate_of_progress_dspecies = multiply({dforward_rate}, forward_reaction_{reaction_index}); // [{reaction_index}][{species_index}] +\n"
                    formatted_text += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)

            for species_index, dbackward_rate in enumerate(backward_rate_derivatives):
                if dbackward_rate == '1':
                    formatted_text += f"        drate_of_progress_dspecies = -divide(forward_reaction_{reaction_index}, equilibrium_constant);\n"
                    formatted_text += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
                elif dbackward_rate != '0':
                    formatted_text += f"        drate_of_progress_dspecies = -multiply({dbackward_rate}, divide(forward_reaction_{reaction_index}, equilibrium_constant));\n"
                    formatted_text += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)

    else:
        if "temperature" in reactions_depend_on[reaction_index]:
            if temperature_jacobian:
                formatted_text += """
        {scalar} drate_of_progress_{reaction_index}_dtemperature =  multiply({forward_rate}, dforward_reaction_{reaction_index}_dtemperature);""".format(reaction_index = reaction_index, 
                    forward_rate = forward_rate,
                    **vars(configuration))
            else:
                formatted_text += """
                // rate_of_progress temperature derivative unused
                """ 
        if ("species" in reactions_depend_on[reaction_index]) or ("pressure" in reactions_depend_on[reaction_index]):
            dforward_rate_dspecies = ''
            for species_index, dforward_rate in enumerate(forward_rate_derivatives):
                if dforward_rate == '1':
                    dforward_rate_dspecies += f"        drate_of_progress_dspecies += forward_reaction_{reaction_index};\n"
                    dforward_rate_dspecies += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
                elif dforward_rate != '0':
                    dforward_rate_dspecies += f"        drate_of_progress_dspecies+= multiply({dforward_rate}, forward_reaction_{reaction_index});\n"
                    dforward_rate_dspecies += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)

            formatted_text += """
        drate_of_progress_dspecies_all_species = 
        scale_gen({forward_rate}, dforward_reaction_{reaction_index}_dspecies);
{all_species}
{dforward_rate_dspecies}
                        """.format(reaction_index = reaction_index, 
                    forward_rate = forward_rate, 
                    dforward_rate_dspecies = dforward_rate_dspecies,
                    all_species = add_to_jacobian_all("drate_of_progress_dspecies_all_species", indexes_of_species_in_reaction, stoichiometric_production),
                    **vars(configuration))
        else:
            formatted_text += """
        //drate_of_progress_dspecies[{reaction_index}] = {{{scalar_cast}(0)}};
""".format(reaction_index = reaction_index, **vars(configuration))
            for species_index, dforward_rate in enumerate(forward_rate_derivatives):
                if dforward_rate == '1':
                    formatted_text += f"        drate_of_progress_dspecies = forward_reaction_{reaction_index};\n"
                    formatted_text += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
                elif dforward_rate != '0':
                    formatted_text += f"        drate_of_progress_dspecies = multiply({dforward_rate}, forward_reaction_{reaction_index});\n"
                    formatted_text += add_to_jacobian("drate_of_progress_dspecies", species_index, indexes_of_species_in_reaction, stoichiometric_production)
    progress_rates_derivatives[reaction_index] = f"        //{reaction}\n"+formatted_text.replace('1 * ','').replace('1.0*','')
def create_equilibrium_constants(stoichiometric_production, reaction_index, indexes_of_species_in_reaction, equilibrium_constants, dequilibrium_constants_dtemperature, configuration, fit_gibbs_reaction = True, temperature_jacobian = False):
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
        if temperature_jacobian:
            dequilibrium_constants_dtemperature[reaction_index] = """
        multiply_chain(exp_gen(-gibbs_reactions[{reaction_index}]), 
                       exp_chain(-gibbs_reactions[{reaction_index}], 
                                 -dgibbs_reactions_dlog_temperature[{reaction_index}]) * dlog_temperature_dtemperature, 
                       {power_term}, 
                       {dpower_term_dtemperature})""".format(power_term = power_term, reaction_index = reaction_index, dpower_term_dtemperature = dpower_term_dtemperature)
        else:
            dequilibrium_constants_dtemperature[reaction_index]="""
        //equilibrium constant temperature derivative is unused
            """
    else:
        print("**WARNING SUPPORT WILL BE DROPPED FOR THIS SOON**")
        equilibrium_constants[reaction_index] = "exp_gen(-({gibbs_sum}) * inv_universal_gas_constant_temperature) * {power_term}".format(gibbs_sum = '+'.join(equilibrium_constant_elements).replace("+-","-"), power_term=power_term)
