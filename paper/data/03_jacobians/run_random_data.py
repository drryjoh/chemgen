#!/usr/bin/env python3
import cantera as ct
import numpy as np
import sys
import chemgen as cg

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 generate_random_states.py <mechanism.yaml> <n_points>")
        sys.exit(1)

    mechanism = sys.argv[1]
    n_points = int(sys.argv[2])

    mech_name = mechanism.split('.')[0]

    random_states = np.load(f"{mech_name}_random_states.npy")
    J = []
    dTdc = []
    for y in random_states:
        J.append(cg.source_jacobian(y[:-1],y[-1]))
        dTdc.append(cg.dtemperature_dspecies(y[:-1],y[-1]))
    
    J = np.array(J)
    dTdc = np.array(dTdc)
    if cg.ignore_temp_dependence():
        np.save(f"J_{mech_name}_itd.npy", J)
    else:
        np.save(f"J_{mech_name}_wtd.npy", J)
        np.save(f"dTdc_{mech_name}.npy", dTdc)


if __name__ == "__main__":
    main()

