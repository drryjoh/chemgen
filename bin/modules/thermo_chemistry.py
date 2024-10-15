import numpy as np
import cantera as ct

pressure_reference = 101325.0

def get_internal_energy_coefficients(gas, order, enthalpy_coefficients):
    internal_internal_energy_coefficients = enthalpy_coefficients.copy()
    internal_internal_energy_coefficients[1,:] = enthalpy_coefficients[1,:] - ct.gas_constant/gas.molecular_weights
    return internal_internal_energy_coefficients

def get_enthalpy_coefficients(gas, order, specific_heat_constant_pressure_species_coefficients):

    enthalpy_coefficients = np.zeros([np.shape(specific_heat_constant_pressure_species_coefficients)[0]+1, np.shape(specific_heat_constant_pressure_species_coefficients)[1]]) 
    for k in range(gas.n_species):
        print(np.double(np.arange(1, np.shape(enthalpy_coefficients)[0])))
        enthalpy_coefficients[1:, k] = specific_heat_constant_pressure_species_coefficients[:, k]/np.double(np.arange(1, np.shape(enthalpy_coefficients)[0]))
    
    Tref = 298.
    gas.TP = Tref, pressure_reference

    T_energy_monomial_sequence = np.power(gas.T, np.linspace(0., order + 1, num = order + 2))
    print(T_energy_monomial_sequence)
    href_polyfit = np.einsum('ij,i->j', enthalpy_coefficients, T_energy_monomial_sequence)
    href_exact = gas.standard_enthalpies_RT * gas.T * ct.gas_constant/gas.molecular_weights # J/kg

    enthalpy_coefficients[0, :] = href_exact - href_polyfit

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

def polyfit_thermodynamics(gas, configuration, order = 4, temperature_min = 200, temperature_max = 8000):
    gas.transport_model = "Mix"
    temperatures = np.linspace(temperature_min, temperature_max, 1000)
    
    #mass specific quantities (units/kg)
    specific_heat_constant_pressure_species_coefficients = get_specific_heat_constant_pressure_species_coefficients(gas, order, temperatures)
    enthalpy_coefficients = get_enthalpy_coefficients(gas, order, specific_heat_constant_pressure_species_coefficients)

    internal_internal_energy_coefficients = get_internal_energy_coefficients(gas, order, enthalpy_coefficients)

    species_specific_heat_text = thermo_fit_text("temperature_monomial_sequence", specific_heat_constant_pressure_species_coefficients, "specific heat", configuration)
    speices_enthalpy_text =      thermo_fit_text("temperature_energy_monomial_sequence", enthalpy_coefficients, "energy", configuration)
    internal_internal_text =      thermo_fit_text("temperature_energy_monomial_sequence", internal_internal_energy_coefficients, "energy", configuration)

    return [["species_specific_heat_constant_pressure_mass_specific",
    "species_internal_energy_mass_specific",
    "species_enthalpy_mass_specific"],
    [species_specific_heat_text,
    speices_enthalpy_text,
    internal_internal_text],
    ["specific_heat", "energy", "energy"]]


def thermo_fit_text(contract_variable, coefficients, thermo_type, configuration, indentation=' '*8):
    thermo_shape = np.shape(coefficients)
    order = thermo_shape[0]
    n_species  = thermo_shape[1]
    content_array = []
    content = f"{indentation}return\n"

    for k in range(n_species):
        coefficients_k = ', '.join(["{scalar_cast}({coefficient})".format(**vars(configuration), coefficient = coefficient) for coefficient in coefficients[:,k]])
        content_array.append("{indentation}contract({contract_variable}, {temperature_monomial_type}{{{coefficients_k}}})".format(indentation = indentation,
        **vars(configuration),
        contract_variable = contract_variable,
        temperature_monomial_type = "{temperature_energy_monomial}".format(**vars(configuration)) if thermo_type == "energy" else "{temperature_monomial}".format(**vars(configuration)),
        coefficients_k = coefficients_k))

    content += '{indentation}{species}{{\n'.format(**vars(configuration), indentation = indentation) + ',\n'.join(content_array)+'};\n'
    return content
