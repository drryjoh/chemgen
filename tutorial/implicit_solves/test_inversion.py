#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import gmres

def source_1pt3(y_n):
    y_n = np.asarray(y_n).flatten()
    s = np.zeros(3)
    s[0] = -0.04 * y_n[0] + 1e4 * y_n[1] * y_n[2]
    s[1] =  0.04 * y_n[0] - 1e4 * y_n[1] * y_n[2] - 3e7 * y_n[1]**2
    s[2] =  3e7 * y_n[1]**2
    return s

def dsource_1pt3_dy(y_n):
    y_n = np.asarray(y_n).flatten()
    s = np.zeros((3,3))
    #s[0] = -0.04 * y_n[0] + 1e4 * y_n[1] * y_n[2]
    s[0,0] = -0.04 
    s[0,1] =  1e4 * y_n[2] 
    s[0,2] =  1e4 * y_n[1]
    #s[1] =  0.04 * y_n[0] - 1e4 * y_n[1] * y_n[2] - 3e7 * y_n[1]**2
    s[1,0] =  0.04
    s[1,1] = -1e4 * y_n[2] - 6e7 * y_n[1]
    s[1,2] = -1e4 * y_n[1]
    #s[2] =  3e7 * y_n[1]**2
    s[2,1] =  6e7 * y_n[1]
    return s

# Initial condition
y0 = np.array([1.0, 0.0, 0.0])
yn = y0.copy()
dt = 1e-3
time_final = 0.3
n_time_steps = int(time_final / dt)
n_newton = 5
ys = [y0]
time = [0.0]
I = np.eye(3)

for t in range(n_time_steps):
    y_guess = yn + 1e-4 * np.ones(3)
    for _ in range(n_newton):
        res = (y_guess - yn) / dt - source_1pt3(y_guess)
        J = I / dt - dsource_1pt3_dy(y_guess)
        dy, info = gmres(J, -res)
        y_guess = y_guess + dy
    yn = y_guess
    ys.append(yn)
    time.append(time[-1] + dt)

ys = np.array(ys)
time = np.array(time)

# Plot y[1] (middle component)
plt.plot(time, ys[:, 1])
plt.xlabel("Time")
plt.ylabel("y[1]")
plt.title("Implicit GMRES Solve of dy/dt = source(y)")
plt.grid(True)
plt.show()
