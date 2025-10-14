#!python3
import numpy as np
import matplotlib.pyplot as plt

# Load stiffness data
stiffness_ffcm2 = np.load("stiffness_FFCM2_model.npy")
stiffness_SanDiego = np.load("stiffness_SanDiego.npy")

# Load temperature data
data = np.loadtxt("SanDiego_data.csv", delimiter=',')
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
ax2.semilogy(time*1e6, stiffness_ffcm2, '-k',markevery=10, label="FFCM2 Stiffness")
ax2.semilogy(time*1e6, stiffness_SanDiego, '-b',markevery=10, label="San Diego Stiffness")

ax2.set_ylabel("Stiffness $||Re(\\lambda)_{max}||/||Re(\\lambda)_{min,\\ne 0}||$", color='b')
ax2.tick_params(axis='y', labelcolor='b')

# Add legends
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

# Show plot
plt.title("Temperature and Stiffness vs. Time")
plt.savefig("c2h4_stiffness.png", dpi=300)

# Create figure and primary axis
fig, ax3 = plt.subplots(figsize=(8, 6))

# Plot Temperature on primary y-axis
ax3.plot(time*1e6, T, '-r', label="Temperature")  
ax3.set_xlabel("Time [$\\mu$s]")
ax3.set_ylabel("Temperature [K]", color='r')
ax3.tick_params(axis='y', labelcolor='r')

# Create secondary y-axis for stiffness values
ax4 = ax3.twinx()  
condition_ffcm2 = np.load("condition_FFCM2_model.npy")
condition_SanDiego = np.load("condition_SanDiego.npy")

ax4.semilogy(time*1e6, condition_ffcm2, '-k',markevery=10, label="FFCM2 Condition #")
ax4.semilogy(time*1e6, condition_SanDiego, '-b',markevery=10, label="San Diego Condition #")


ax4.set_ylabel("Condition $||Re(\\lambda)_{max}||/||Re(\\lambda)_{min,\\ne 0}||$", color='b')
ax4.tick_params(axis='y', labelcolor='b')

# Add legends
ax3.legend(loc="upper left")
ax4.legend(loc="upper right")

# Show plot
plt.title("Temperature and Condition vs. Time")
plt.savefig("c2h4_condition.png", dpi=300)
plt.show()

