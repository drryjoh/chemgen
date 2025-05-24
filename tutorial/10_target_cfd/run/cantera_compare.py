import cantera as ct
import numpy as np
gas = ct.Solution("chemkin/chem.yaml")

X_string = "H2:1, N2:3.76, O2:1"
p = 202650
T = 1000

gas.TPX = T, p, X_string 
npr = gas.net_production_rates * gas.molecular_weights
rhoU  = gas.int_energy_mass * gas.density
print(rhoU)
print(gas.concentrations)