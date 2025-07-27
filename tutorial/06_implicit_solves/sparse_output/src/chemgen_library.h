#include "types_inl.h"

int species_index_gen(const char* name);


double
density_from_massfractions_pressure_temperature(const Species& massfractions,
                                                const double& pressure,
                                                const double& temperature) ;


Species
scale_gen(const double& a, const Species& b_s) ;


ChemicalState
scale_gen(const double& a, const ChemicalState& b_s) ;


Species
operator+(const Species& a_s, const Species& b_s) ;


SpeciesJacobian
operator+(const SpeciesJacobian& a, const SpeciesJacobian& b) ;


Species
operator-(const Species& a_s, const Species& b_s_positive) ;


ChemicalState
operator+(const ChemicalState& a_s, const ChemicalState& b_s) ;


ChemicalState
operator-(const ChemicalState& a_s, const ChemicalState& b_s_positive) ;


Species
operator*(const double& s, const Species& a) ;


ChemicalState
operator*(const double& s, const ChemicalState& a) ;


Species
operator*(const Species& a_s, const Species& b_s) ;


double
norm2(const Species& x) ;


double
divide(const double& a, const double& b) ;


Species
divide(const Species& a, const Species& b) ;


double
sqrt_gen(const double& a) ;


double
scale_gen(const double& a, const double& b) ;

Species inv_molecular_weights_();


Species
concentrations_from_molefractions_pressure_temperature(const Species& molefractions,
                                                       const double& pressure,
                                                       const double& temperature) ;


double
internal_energy_volume_specific(const Species& species,
                                const double& temperature) ;


ChemicalState
set_chemical_state(const double& energy,
                   const Species& species) ;


double
temperature(const ChemicalState& chemical_state) ;


double
temperature(double energy, const ChemicalState& chemical_state) ;


ChemicalState
source(const ChemicalState& chemical_state, double temperature) ;


ChemicalState
source(const ChemicalState& chemical_state) ;


SparseMatrix<double>
source_jacobian(const ChemicalState& chemical_state, double temperature, const double& scaling_factor=1, const double& diagonal_add=0) ;


Species
get_species(const ChemicalState& y) ;


double error_norm(const ChemicalState& x,
                             const ChemicalState& x_ref,
                             double abs_tolerance,
                             double rel_tolerance) ;