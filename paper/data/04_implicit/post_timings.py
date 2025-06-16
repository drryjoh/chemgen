#!python3
import numpy as np
'''
[Backward Euler] Time elapsed: 0.225836 seconds
[SDIRK2] Time elapsed: 0.385753 seconds
[ROSENBROC] Time elapsed: 0.0907028 seconds
[YASS] Time elapsed: 0.0585083 seconds
[RK4] Time elapsed: 0.326525 seconds
[SDIRK4] Time elapsed: 0.261711 seconds
'''
time = np.array([0.222798, 0.385753, 0.0907028, 0.0585083, 0.309098 , 0.261711])
cantera_time =  0.0586651
time_norm = time/cantera_time
imp_types = "Backward Euler, SDIRK2, ROSENBROC, YASS, RK4, SDIRK4"
imp_types_array = imp_types.split(",")
for t, it in zip(time_norm, imp_types_array):
    print(f"{it} : {t:3.2f}")
