#!/usr/bin/env python3
import cantera as ct
import numpy as np
import sys

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
    print(f"mech = {mechanism}")

    gas = ct.Solution(mechanism)
    random_states = []

    for _ in range(n_points):
        random_states.append(get_random_TPX(gas))

    random_states = np.array(random_states)
    mech_name = mechanism.split('.')[0]

    np.save(f"{mech_name}_random_states.npy", random_states)
    print(f"Saved {n_points} random states to {mech_name}_random_states.npy")

if __name__ == "__main__":
    main()

