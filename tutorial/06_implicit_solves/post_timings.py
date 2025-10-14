#!python3
import numpy as np
'''
[Backward Euler] Time elapsed: 0.191271 seconds
[SDIRK2] Time elapsed: 0.0758483 seconds
[ROSENBROC] Time elapsed: 0.106801 seconds
[YASS] Time elapsed: 0.075036 seconds
[RK4] Time elapsed: 0.286059 seconds
[SDIRK4] Time elapsed: 0.0687864 seconds
'''
time = np.array([0.191271, 0.0758483, 0.106801, 0.075036, 0.309098 , 0.0687864])
cantera_time =   0.241663
time_norm = cantera_time/time
imp_types = "Backward Euler, SDIRK2, ROSENBROC, YASS, RK4, SDIRK4"
imp_types_array = imp_types.split(",")
for t, it in zip(time_norm, imp_types_array):
    print(f"{it} : {t:3.2f}")
