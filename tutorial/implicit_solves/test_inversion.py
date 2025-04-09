#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

import time
from implicit_schemes import backwards_euler, sdirk2, sdirk4, rosenbrock2, seulex_adaptive

def timeit_solver(name, func, *args, **kwargs):
    print(f"Running {name}...")
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    print(f"{name} completed in {elapsed:.4f} seconds\n")
    return result, elapsed

# Initial condition
y0 = np.array([1.0, 0.0, 0.0])
dt = 1e-3
time_final = 0.3
n_time_steps = int(time_final / dt)

(time_be, y_be), t_be = timeit_solver("Backward Euler", backwards_euler, y0, dt, n_time_steps, n_newton=5)
(time_sdirk2, y_sdirk2), t_sdirk2 = timeit_solver("SDIRK2", sdirk2, y0, dt, n_time_steps, n_newton=5)
(time_sdirk4, y_sdirk4), t_sdirk4 = timeit_solver("SDIRK4", sdirk4, y0, dt, n_time_steps, n_newton=5)
(time_be_h2, y_be_h2), t_be_h2 = timeit_solver("Backward Euler (h/2)", backwards_euler, y0, dt/2.0, 2*n_time_steps, n_newton=5)
(time_sdirk2_h2, y_sdirk2_h2), t_sdirk2_h2 = timeit_solver("SDIRK2 (h/2)", sdirk2, y0, dt/2, 2*n_time_steps, n_newton=5)
(time_ros2, y_ros2), t_ros2 = timeit_solver("Rosenbrock2", rosenbrock2, y0, dt, n_time_steps)
(time_ros2_h2, y_ros2_h2), t_ros2_h2 = timeit_solver("Rosenbrock2", rosenbrock2, y0, dt/2, 2*n_time_steps)
(time_sdirk4_h10, y_sdirk4_h10), t_sdirk4_h10 = timeit_solver("SDIRK4", sdirk4, y0, dt/10, n_time_steps*10, n_newton=5)
(time_seulex, y_seulex), t_seulex = timeit_solver("SEULEX", seulex_adaptive, y0, dt, time_final)

plt.plot(time_sdirk4_h10, y_sdirk4_h10[:, 1], '-k', label=f'SDIRK4, dt/10\ntime = {t_sdirk4_h10:.2e}', mfc = "white", lw=4)
plt.plot(time_be, y_be[:, 1], '-or', label=f'Backward Euler\ntime = {t_be:.2e}', mfc = "white")
plt.plot(time_sdirk2, y_sdirk2[:, 1], '-ok', label=f'SDIRK2\ntime = {t_sdirk2:.2e}', mfc = "white")
plt.plot(time_ros2, y_ros2[:, 1], '-og', label=f'Rosenbrock2\ntime = {t_ros2:.2e}', mfc = "white")
plt.plot(time_sdirk4, y_sdirk4[:, 1], '-ob', label=f'SDIRK4\ntime = {t_sdirk4:.2e}', mfc = "white")
plt.plot(time_be_h2, y_be_h2[:, 1], '-^r', label=f'Backward Euler, dt/2\ntime = {t_be_h2:.2e}', mfc = "white")
plt.plot(time_sdirk2_h2, y_sdirk2_h2[:, 1], '-^k', label=f'SDIRK2, dt/2\ntime = {t_sdirk2_h2:.2e}', mfc = "white")
plt.plot(time_ros2_h2, y_ros2_h2[:, 1], '-^g', label=f'Rosenbrock2, dt/2\ntime = {t_ros2_h2:.2e}', mfc = "white")
plt.plot(time_seulex, y_seulex[:, 1], '-^r', label=f'SEULEX, dt/2\ntime = {t_seulex:.2e}', mfc = "white")

np.save("time_sdirk4_h10.npy", time_sdirk4_h10)
np.save("y_sdirk4_h10.npy", y_sdirk4_h10)

np.save("time_be.npy", time_be)
np.save("y_be.npy", y_be)

np.save("time_sdirk2.npy", time_sdirk2)
np.save("y_sdirk2.npy", y_sdirk2)

np.save("time_ros2.npy", time_ros2)
np.save("y_ros2.npy", y_ros2)

np.save("time_sdirk4.npy", time_sdirk4)
np.save("y_sdirk4.npy", y_sdirk4)

np.save("time_be_h2.npy", time_be_h2)
np.save("y_be_h2.npy", y_be_h2)

np.save("time_sdirk2_h2.npy", time_sdirk2_h2)
np.save("y_sdirk2_h2.npy", y_sdirk2_h2)

np.save("time_ros2_h2.npy", time_ros2_h2)
np.save("y_ros2_h2.npy", y_ros2_h2)

plt.ylim([3.636e-5, 3.65e-5])
plt.xlim([0.0025, 0.02])
plt.legend(loc=1, ncol=2, fontsize=8)
plt.xlabel("Time")
plt.ylabel("$y_2$")
plt.title("Implicit Solve of dy/dt = S(y)")
plt.grid(True)
plt.show()
