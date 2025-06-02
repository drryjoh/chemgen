
const int n_species = 96;
const int n_reactions = 1054;
const int n_order_thermo = 7 + 1;
const int n_chemical_state = 96 + 1;
// Using alias for the array type (for example, an array of double values)
using Species = std::array<double, n_species>;
using Reactions = std::array<double, n_reactions>;
using TemperatureMonomial = std::array<double, n_order_thermo>;
using TemperatureEnergyMonomial = std::array<double, n_order_thermo + 1>;
using TemperatureGibbsMonomial = std::array<double, n_order_thermo + 2>;
using ThermoTable = std::array<TemperatureEnergyMonomial, n_species>;
using ChemicalState = std::array<double, n_species + 1>;
using SpeciesJacobian = std::array<Species, n_species>;
using ReactionSpecies = std::array<std::array<double, n_species>, n_reactions>;

