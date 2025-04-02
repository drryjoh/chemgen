import numpy as np
import matplotlib.pyplot as plt
time_py = np.load("time_sdirk4_h10.npy")
sdirk_py = np.load("y_sdirk4_h10.npy")
data = np.loadtxt("backward_euler.csv", delimiter=",", skiprows=1)

# Extract columns
time = data[:, 0]
y0 = data[:, 1]
y1 = data[:, 2]
y2 = data[:, 3]

# Plot
plt.plot(time, y1, label='y1 from C++ Backward Euler')
plt.plot(time_py, sdirk_py[:,1], label='y1 from python Backward Euler')
plt.xlabel('Time')
plt.ylabel('y1')
plt.legend()
plt.grid(True)
plt.title("C++ Backward Euler Result")
plt.show()
