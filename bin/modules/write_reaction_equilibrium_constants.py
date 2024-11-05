def write_equilibrium_constants(file, equilibrium_constants, configuration):
    equilibrium_constant_evaluations = ''
    indentation=' '*8
    for i, equilibrium_constant_i in enumerate(equilibrium_constants):
        equilibrium_constant_evaluations += f"{indentation}equilibrium_constants[{i}] = {equilibrium_constant_i};\n"
    file.write("""
    {device_option}
    {reactions_function} equilibrium_constants({scalar_parameter} temperature) {const_option} 
    {{
        {species} gibbs_free_energies = species_gibbs_energy_mole_specific(temperature);
        {reactions} equilibrium_constants = {{}};
        {scalar} inv_universal_gas_constant_temperature  = inv(universal_gas_constant() * temperature);
{equilibrium_constant_evaluations}
        return equilibrium_constants;
    }}
""".format(**vars(configuration),
           equilibrium_constant_evaluations = equilibrium_constant_evaluations))