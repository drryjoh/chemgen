import numpy as np
import matplotlib.pyplot as plt

# Load stiffness data
stiffness_ffcm2 = np.load("stiffness_ffcm2_h2.npy")
stiffness_burke = np.load("stiffness_burke.npy")
stiffness_h2o2 = np.load("stiffness_h2o2.npy")

# Load temperature data
data = np.loadtxt("burke_data.csv", delimiter=',')
time = data[:, 0]  # Time
T = data[:, 1]  # Temperature

# Create figure and primary axis
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot Temperature on primary y-axis
ax1.plot(time*1e6, T, '-r', label="Temperature")  
ax1.set_xlabel("Time [$\\mu$s]")
ax1.set_ylabel("Temperature [K]", color='r')
ax1.tick_params(axis='y', labelcolor='r')

# Create secondary y-axis for stiffness values
ax2 = ax1.twinx()  
ax2.semilogy(time*1e6, stiffness_ffcm2, '-b', label="FFCM2")
ax2.semilogy(time*1e6, stiffness_burke, '--b', label="Burke")
ax2.semilogy(time*1e6, stiffness_h2o2, '-.b', label="GRI 3.0")

ax2.set_ylabel("Stiffness", color='b')
ax2.tick_params(axis='y', labelcolor='b')

# Add legends
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

# Show plot
plt.title("Temperature and Stiffness vs. Time")
plt.show()

