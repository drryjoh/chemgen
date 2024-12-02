import numpy as np
import cantera as ct
reactions = [ 184,  538,  584,  770,  911,  945,  990, 1008]

gas = ct.Solution("FFCM2_model.yaml")

for reaction_index in reactions:
    print(gas.reactions()[reaction_index].reaction_type)
