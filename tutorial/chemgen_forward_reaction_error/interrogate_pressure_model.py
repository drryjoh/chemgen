import numpy as np
import cantera as ct
gas = ct.Solution("test.yaml")
for reaction in gas.reactions():
    for p, rate  in reaction.rate.rates:
        print(p)
