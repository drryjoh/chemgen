import numpy as np
import cantera as ct


pressure_reference = 101325.0
def get_gibbs_energy_coefficients(gas, order, specific_heat_coefficients, enthalpy_coefficients, entropy_coefficients):
    specific_heat_shape = np.shape(specific_heat_coefficients)
    gibbs_energy_coefficients = np.zeros([specific_heat_shape[0] + 2, specific_heat_shape[1]])

    h_ref = enthalpy_coefficients[0,:]
    s_ref = entropy_coefficients[0, :]

    gibbs_energy_coefficients[0, :] = h_ref
    gibbs_energy_coefficients[1, :] =  specific_heat_coefficients[0,:] - s_ref
    for k in range(2, order+1):
        gibbs_energy_coefficients[k, :] = (1/k - 1/(k-1)) * specific_heat_coefficients[k-1,:]

    gibbs_energy_coefficients[-1, :] = -specific_heat_coefficients[0,:] #TlnT term

    for k in range(gas.n_species):
        gibbs_energy_coefficients[:,k]  = gibbs_energy_coefficients[:,k] * gas.molecular_weights[k]


    return gibbs_energy_coefficients


def get_entropy_coefficients(gas, order, internal_energy_coefficients, specific_heat_constant_pressure_species_coefficients):
    entropy_coefficients = np.zeros_like(internal_energy_coefficients)
            
    entropy_log_coefficient = specific_heat_constant_pressure_species_coefficients[0,:] # a_0
    for k in range(gas.n_species):
        for i, a in enumerate(specific_heat_constant_pressure_species_coefficients[:,k]):
            if i==0:
                pass
            else:
                entropy_coefficients[i,k] = a/(i)
    entropy_coefficients[-1,:] = entropy_log_coefficient
    #assure that 298 matches
    Tref = 298.
    gas.TP = Tref, pressure_reference
    T_entropy_monomial_sequence = np.power(gas.T, np.linspace(0., order+1, num = order + 2))
    T_entropy_monomial_sequence[-1] = np.log(Tref)


    sref_polyfit = []
    for k in range(gas.n_species):
        sref_polyfit.append(np.sum(entropy_coefficients[:,k] * T_entropy_monomial_sequence))
    sref_polyfit = np.array(sref_polyfit)
    sref_exact = gas.standard_entropies_R * ct.gas_constant/gas.molecular_weights # J/kg
    entropy_coefficients[0, :] += sref_exact - sref_polyfit

    return entropy_coefficients

def get_internal_energy_coefficients(gas, order, enthalpy_coefficients):
    internal_internal_energy_coefficients = enthalpy_coefficients.copy()
    internal_internal_energy_coefficients[1,:] = enthalpy_coefficients[1,:] - ct.gas_constant/gas.molecular_weights
    return internal_internal_energy_coefficients

def check_against_temperature(gas, temperature_min, temperature_max, n_samples, enthalpy_coefficients, order):
    import matplotlib.pyplot as plt
    temperatures = np.linspace(temperature_min, temperature_max, n_samples)
    h_exact = []
    h_refit = []
    h_exact = np.zeros([gas.n_species, n_samples])
    h_fit = np.zeros([gas.n_species, n_samples])
    for i, temperature in enumerate(temperatures):
        gas.TP = temperature, pressure_reference
        T_energy_monomial_sequence = np.power(gas.T, np.linspace(0., order + 1, num = order + 2))
        href_exact = gas.standard_enthalpies_RT * gas.T * ct.gas_constant/gas.molecular_weights # J/kg
        for k in range(gas.n_species):
            h_exact[k,i] = href_exact[k]
            h_fit[k,i] = np.sum(enthalpy_coefficients[:,k] * T_energy_monomial_sequence)
    for k in range(gas.n_species):
        plt.figure()
        plt.plot(temperatures, h_exact[k,:],'-r')
        plt.plot(temperatures, h_fit[k,:],'--k')
    plt.show()

def get_enthalpy_coefficients(gas, order, specific_heat_coefficients):
    specific_heat_shape = np.shape(specific_heat_coefficients)
    enthalpy_coefficients = np.zeros([specific_heat_shape[0]+1, specific_heat_shape[1]]) 

    for i in range(1,order+1):
        enthalpy_coefficients[i,:] = specific_heat_coefficients[i-1,:]/(i)

    #assure that 298 matches
    Tref = 298.
    gas.TP = Tref, pressure_reference
    T_energy_monomial_sequence = np.power(gas.T, np.linspace(0., order + 1, num = order + 2))

    href_polyfit = []
    for k in range(gas.n_species):
        href_polyfit.append(np.sum(enthalpy_coefficients[:,k] * T_energy_monomial_sequence))
    href_polyfit = np.array(href_polyfit)
    href_exact = gas.standard_enthalpies_RT * gas.T * ct.gas_constant/gas.molecular_weights # J/kg
    enthalpy_coefficients[0, :] = href_exact - href_polyfit
    
    #check_against_temperature(gas, 200, 5000, 100, enthalpy_coefficients, order)
    return enthalpy_coefficients

def get_specific_heat_constant_pressure_species_coefficients(gas, order, temperatures):
    specific_heats_constant_pressure = []

    for temperature in temperatures:
        gas.TP = temperature, pressure_reference
        specific_heats_constant_pressure.append(gas.standard_cp_R)
    
    specific_heats_constant_pressure = np.array(specific_heats_constant_pressure)

    #polyfit
    specific_heat_constant_pressure_species_coefficients = np.polyfit(temperatures, 
    specific_heats_constant_pressure * ct.gas_constant / gas.molecular_weights,
    order) # mass-specific, dimensional

    polynomial = np.poly1d(specific_heat_constant_pressure_species_coefficients[:,0])

    #reverse order
    for k in range(gas.n_species):
        specific_heat_constant_pressure_species_coefficients[:,k] = specific_heat_constant_pressure_species_coefficients[::-1,k]

    return specific_heat_constant_pressure_species_coefficients

def get_gibbs_reaction_coefficients(gas, order, temperatures):
    gibbs_reaction = []

    for temperature in temperatures:
        gas.TP = temperature, pressure_reference
        gibbs_reaction.append(gas.delta_standard_gibbs/ct.gas_constant/temperature)
    
    gibbs_reaction = np.array(gibbs_reaction)

    #polyfit
    gibbs_reaction_coefficients = np.polyfit(np.log(temperatures), gibbs_reaction, order)

    #reverse order
    for k in range(gas.n_reactions):
        gibbs_reaction_coefficients[:,k] = gibbs_reaction_coefficients[::-1,k]

    return gibbs_reaction_coefficients

def check_gibbs_against_temperature(gas, temperature_min, temperature_max, n_samples, gibbs_coefficients, order):
    import matplotlib.pyplot as plt
    temperatures = np.linspace(temperature_min, temperature_max, n_samples)
    gibbs_exact = np.zeros([gas.n_reactions, n_samples])
    gibbs_fit = np.zeros([gas.n_reactions, n_samples])

    log_T_energy_monomial_sequence = np.power(np.log(1800), np.linspace(0., order, num = order + 1))
    print(log_T_energy_monomial_sequence)
    print(gibbs_coefficients)

    for i, temperature in enumerate(temperatures):
        gas.TP = temperature, pressure_reference
        log_T_energy_monomial_sequence = np.power(np.log(gas.T), np.linspace(0., order, num = order + 1))
        
        gibbs_ct = gas.delta_standard_gibbs/ct.gas_constant/temperature
        for k in range(gas.n_reactions):
            gibbs_exact[k,i] = gibbs_ct[k]
            gibbs_fit[k,i] = np.sum(gibbs_coefficients[:,k] * log_T_energy_monomial_sequence)
    for k in range(gas.n_reactions):
        plt.figure()
        plt.plot(temperatures, gibbs_exact[k,:],'-r')
        plt.plot(temperatures, gibbs_fit[k,:],'--k')
    plt.show()
def polyfit_thermodynamics(gas, configuration, order = 4, temperature_min = 200, temperature_max = 8000, transport = False):
    if transport:
        gas.transport_model = "Mix"
    temperatures = np.linspace(temperature_min, temperature_max, 1000)
    
    #mass specific quantities (units/kg)
    specific_heat_constant_pressure_species_coefficients = get_specific_heat_constant_pressure_species_coefficients(gas, order, temperatures)
    enthalpy_coefficients = get_enthalpy_coefficients(gas, order, specific_heat_constant_pressure_species_coefficients)
    internal_energy_coefficients = get_internal_energy_coefficients(gas, order, enthalpy_coefficients)
    species_entropy_coefficients = get_entropy_coefficients(gas, order, internal_energy_coefficients, specific_heat_constant_pressure_species_coefficients)
    gibbs_energy_coefficients = get_gibbs_energy_coefficients(gas, order, specific_heat_constant_pressure_species_coefficients, enthalpy_coefficients, species_entropy_coefficients)
    gibbs_reaction_energy_coefficients = get_gibbs_reaction_coefficients(gas, order, temperatures)

    species_specific_heat_text = thermo_fit_text("temperature_monomial_sequence", specific_heat_constant_pressure_species_coefficients, "default", configuration)
    species_enthalpy_text  =      thermo_fit_text("temperature_energy_monomial_sequence", enthalpy_coefficients, "energy", configuration)
    internal_internal_text =      thermo_fit_text("temperature_energy_monomial_sequence", internal_energy_coefficients, "energy", configuration)
    species_entropy_text   =      thermo_fit_text("temperature_entropy_monomial_sequence", species_entropy_coefficients, "energy", configuration)
    gibbs_energy_text   =      thermo_fit_text("temperature_gibbs_monomial_sequence", gibbs_energy_coefficients, "gibbs", configuration)
    gibbs_energy_reaction_text   =      thermo_fit_text("log_temperature_monomial_sequence", gibbs_reaction_energy_coefficients, "default", configuration, return_type = "{reactions}")

    return [["species_specific_heat_constant_pressure_mass_specific",
    "species_enthalpy_mass_specific",
    "species_internal_energy_mass_specific",
    "species_entropy_mass_specific",
    "species_gibbs_energy_mole_specific",
    "gibbs_reaction"],
    [species_specific_heat_text,
    species_enthalpy_text,
    internal_internal_text,
    species_entropy_text,
    gibbs_energy_text,
    gibbs_energy_reaction_text],
    ["specific_heat", "energy", "energy", "entropy", "gibbs", "gibbs_reaction"]]


def thermo_fit_text(contract_variable, coefficients, thermo_type, configuration, indentation=' '*8, return_type = '{species}'):
    thermo_shape = np.shape(coefficients)
    order = thermo_shape[0]
    n_species  = thermo_shape[1]
    content_array = []
    content = f"{indentation}return\n"
    return_type = return_type.format(**vars(configuration))

    for k in range(n_species):
        coefficients_k = ', '.join(["{scalar_cast}({coefficient})".format(**vars(configuration), coefficient = coefficient) for coefficient in coefficients[:,k]])
        content_array.append("{indentation}contract({contract_variable}, {temperature_monomial_type}{{{coefficients_k}}})".format(indentation = indentation,
        **vars(configuration),
        contract_variable = contract_variable,
        temperature_monomial_type = "{temperature_monomial}".format(**vars(configuration)) if thermo_type == "default" else f"{{temperature_{thermo_type}_monomial}}".format(**vars(configuration)) , 
        coefficients_k = coefficients_k))

    content += f'{indentation}{return_type}'+'{{\n'.format(**vars(configuration), indentation = indentation) + ',\n'.join(content_array)+'};\n'
    return content
