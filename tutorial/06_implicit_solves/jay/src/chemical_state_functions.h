#define CHEMGEN_INTERNAL_ENERGY_EQUATION


double
temperature(const ChemicalState& chemical_state) 
{
    return
#if defined(CHEMGEN_INTERNAL_ENERGY_EQUATION)
    temperature(get_energy(chemical_state), get_species(chemical_state));
#elif defined(CHEMGEN_TEMPERATURE_EQUATION)
    get_energy(chemical_state);
#else
#error
#endif
}


double
temperature(double energy, const ChemicalState& chemical_state) 
{
    return
#if defined(CHEMGEN_INTERNAL_ENERGY_EQUATION)
    temperature(energy, get_species(chemical_state));
#elif defined(CHEMGEN_TEMPERATURE_EQUATION)
    energy;
#else
#error
#endif
}


ChemicalState 
source(const ChemicalState& chemical_state, double temperature) 
{
    return 
    source(get_species(chemical_state), temperature);
}


ChemicalState
source(const ChemicalState& chemical_state) 
{
    return
    source(get_species(chemical_state), temperature(chemical_state));
}


SparseMatrix<double>
source_jacobian(const ChemicalState& chemical_state, double temperature, const double& scaling_factor=1, const double& diagonal_add=0) 
{
    return
    source_jacobian(get_species(chemical_state), temperature, scaling_factor, diagonal_add);
}

#undef CHEMGEN_INTERNAL_ENERGY_EQUATION