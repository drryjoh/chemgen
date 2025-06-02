
double 
temperature(const ChemicalState& chemical_state)  
{
    return temperature(get_energy(chemical_state), get_species(chemical_state));
}


ChemicalState 
source(const ChemicalState& chemical_state)  
{
    return 
    set_chemical_state(double(0),
                       source(get_species(chemical_state), temperature(chemical_state)));
}