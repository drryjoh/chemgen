import numpy as np
import cantera as ct

pressure_reference = 101325.0

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
    print("**test**")
    print(polynomial(1800))
    print(specific_heat_constant_pressure_species_coefficients[:,0])
    for k in range(gas.n_species):
        specific_heat_constant_pressure_species_coefficients[:,k] = specific_heat_constant_pressure_species_coefficients[::-1,k]

    return specific_heat_constant_pressure_species_coefficients

def polyfit_thermodynamics(gas, order = 4, temperature_min = 200, temperature_max = 8000):
    gas.transport_model = "Mix"
    temperatures = np.linspace(temperature_min, temperature_max, 1000)
    return get_specific_heat_constant_pressure_species_coefficients(gas, order, temperatures)

def thermo_fit_text(coefficients, configuration, indentation=' '*8):
    thermo_shape = np.shape(coefficients)
    order = thermo_shape[0]
    n_species  = thermo_shape[1]
    content_array = []
    content = f"{indentation}return\n"
    for k in range(n_species):
        coefficients_k = ', '.join(["{scalar_cast}({coefficient})".format(**vars(configuration), coefficient = coefficient) for coefficient in coefficients[:,k]])
        content_array.append("{indentation}contract(temperature_monomial_sequence, {temperature_monomial}{{{coefficients_k}}})".format(indentation = indentation,
        **vars(configuration),
        coefficients_k = coefficients_k))

    content += '{indentation}{species}{{\n'.format(**vars(configuration), indentation = indentation) + ',\n'.join(content_array)+'};\n'
    return content
