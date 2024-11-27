import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
data = np.loadtxt("l2_norm_results.csv", delimiter=',', skiprows = 1)
data_no_temp = data[1:]
reaction_data = data_no_temp[:-1]
chemgen = reaction_data[::2]
cantera = reaction_data[1::2]
difference = [(cg - ct) / ct if ct != 0 else 0 for cg, ct in zip(chemgen, cantera)]
gas = ct.Solution("double_a.yaml")

gas.TPX = 2000, 101325.0, "OH:0.05, H:0.05, H2O:0.9"#get_random_TPX(gas)
print(gas.species_names)
print(gas.X)
print(gas.concentrations)
print(gas.T)
print(gas.P)
print(gas.forward_rate_constants)
k_f = gas.forward_rate_constants[0]
print(gas.equilibrium_constants)
k_eq = gas.equilibrium_constants[0]
print(f"reverse_rate_constants = {gas.reverse_rate_constants[0]}")
k_r= gas.reverse_rate_constants[0]
print(f"k_f/k_eq/k_r = {k_f/k_eq/k_r}")
C = gas.concentrations
omega = gas.net_production_rates

production_rate_a = k_f*C[-1]*C[-1] + k_r * C[0] * C[1] * C[-1]
production_rate_b = k_f*C[-1] + k_r * C[0] * C[1]
print(production_rate_a)
print(production_rate_b)
print(f"cantera evaluation:  {omega}")
print(f"reaction a approach: {np.array([-1,1,1]) * production_rate_a}")
print(f"reaction b approach: {np.array([-1,1,1]) * production_rate_b}")

print(gas.product_stoich_coeff("H2O",0))
print(gas.reactant_stoich_coeff("H2O",0))
print(gas.reactant_stoich_coeffs)


for k, diff in enumerate(difference):
    if diff>1:
        print(f"Something is up with reaction {k+1}: {diff} {chemgen[k]} {cantera[k]}")
        print(gas.reactions()[k])
        
plt.show()
