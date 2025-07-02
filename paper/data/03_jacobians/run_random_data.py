#!/usr/bin/env python3
import cantera as ct
import numpy as np
import sys
import chemgen as cg
def get_random_TPX(gas):
    species_array = np.random.uniform(0, 1, len(gas.species_names))
    species_array /= species_array.sum()

    T = 1000 + 1500 * np.random.random()
    p = 10132.5 + 101325.0 * 10.9 * np.random.random()

    C = species_array * p / (ct.gas_constant * T)
    return np.hstack((C, T))  # return a flat array with species and temperature

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 generate_random_states.py <mechanism.yaml> <n_points>")
        sys.exit(1)

    mechanism = sys.argv[1]
    n_points = int(sys.argv[2])

    mech_name = mechanism.split('.')[0]

    random_states = np.load(f"{mech_name}_random_states.npy")
    J = []
    for y in random_states:
        J.append(cg.source_jacobian(y[:-1],y[-1]))
    
    J = np.array(J)
    if cg.ignore_temp_dependence():
        np.save(f"J_{mech_name}_itd.npy", J)
    else:
        np.save(f"J_{mech_name}_wtd.npy", J)

if __name__ == "__main__":
    main()

