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

def backwards_euler(y0, dt, n_time_steps, n_newton = 5):
    yn = y0.copy()
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
    return time, ys

y0 = np.array([1.0, 0.0, 0.0])
dt = 1e-2
time_final = 0.3
n_time_steps = int(time_final / dt)

time_be, y_be = backwards_euler(y0, dt, n_time_steps, n_newton = 5)
time_be_h2, y_be_h2 = backwards_euler(y0, dt/2.0, 2 * n_time_steps, n_newton = 5)
time_be_h4, y_be_h4 = backwards_euler(y0, dt/4.0, 4 * n_time_steps, n_newton = 5)
time_be_h8, y_be_h8 = backwards_euler(y0, dt/8.0, 8 * n_time_steps, n_newton = 5)
# Plot y[1] (middle component)
plt.plot(time_be, y_be[:, 1])
plt.plot(time_be_h2, y_be_h2[:, 1])
plt.plot(time_be_h4, y_be_h4[:, 1])
plt.plot(time_be_h8, y_be_h8[:, 1])
plt.ylim([3.5e-5, 3.7e-5])
plt.xlabel("Time")
plt.ylabel("y[1]")
plt.title("Implicit GMRES Solve of dy/dt = source(y)")
plt.grid(True)
plt.show()
