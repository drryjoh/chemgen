#!python3
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Load the reaction mechanism
# gas = ct.Solution("mechanisms/ic8-874-6864.yaml")
gas = ct.Solution("mechanisms/mmc5-7171-38324.yaml")

# Define initial conditions
test_conditions = {
    "temperature": 1800,  # K
    "pressure": 101325.0,  # Pa
    "species": 
    {
        # "O2": 0.14788732394366194,
        # "N2": 0.5563380281690141,
        # "H2": 0.2957746478873239

        # "O2": 0.20736086175942545,
        # "N2": 0.7800718132854579,
        # "POSF10325": 0.012567324955116694

        # "O2": 0.20653029110936266,
        # "N2": 0.7769472856018884,
        # "IC8H18": 0.016522423288749016

        # "O2": 0.2060660124888492,
        # "N2": 0.7752007136485282,
        # "NC7H16": 0.01873327386262266

        "O2": 0.20856398567241938,
        "N2": 0.7845978508629112,
        "c20h42-2": 0.006838163464669487
    }
}

gas.TPX = (
    test_conditions["temperature"],
    test_conditions["pressure"],
    test_conditions["species"]
)

# Create a reactor and insert the gas
reactor = ct.IdealGasReactor(gas)
network = ct.ReactorNet([reactor])

# Define simulation time (in seconds)
time_end = 2e-5
dt_small = 1e-7
n_steps = int(time_end/dt_small)
#time_end = 200 * 2e-7  # Convert ns to seconds
#n_steps = 200  # Number of time steps
time = np.linspace(0, time_end, n_steps)

temperature = []
data = []
import time as clock 
start_time = clock.time()
for t in time:
    network.advance(t)
    temperature.append(reactor.T)
    data.append([t, reactor.T])
end_time = clock.time()
print(f"Elapsed time: {end_time - start_time:.6f} seconds")
data = np.array(data)
# Save results to file

import warnings
warnings.filterwarnings('ignore')

# Plot results
plt.plot(data[:, 0]*1000.0, data[:, 1],'-r', label = "Cantera", lw=4)
d = np.loadtxt("backward_euler.txt")
try:
    plt.plot(d[:, 0]*1000.0, d[:, 1],'-ok', label = "ChemGen Backward Euler", markevery=int(len(d[:, 0])/10))
except IndexError:
    pass
d = np.loadtxt("sdirk2.txt")
try:
    plt.plot(d[:, 0]*1000.0, d[:, 1],'-ob', label = "ChemGen SDIRK2", markevery=int(len(d[:, 0])/10))
except IndexError:
    pass
d = np.loadtxt("sdirk4.txt")
try:
    plt.plot(d[:, 0]*1000.0, d[:, 1],'-og', label = "ChemGen SDIRK4", markevery=int(len(d[:, 0])/10))
except IndexError:
    pass
d = np.loadtxt("ros.txt")
try:
    plt.plot(d[:, 0]*1000.0, d[:, 1],'-o',color='purple', label = "ChemGen Rosenbroc", markevery=int(len(d[:, 0])/10))
except IndexError:
    pass
d = np.loadtxt("yass.txt")
try:
    plt.plot(d[:, 0]*1000.0, d[:, 1],'-d',color='orange', label = "ChemGen YASS", markevery=int(len(d[:, 0])/10))
except IndexError:
    pass
d = np.loadtxt("rk4.txt")
try:
    plt.plot(d[:, 0]*1000.0, d[:, 1],'--ok', label = "ChemGen RK4",markevery=int(len(d[:, 0])/10))
except IndexError:
    pass

plt.legend()
plt.xlabel("Time ($\\mu$s)")
plt.ylabel("Temperature (K)")
# plt.ylim([300,4000])
plt.title("Temperature Evolution in Homogeneous Reactor")
plt.savefig("rk4.png",dpi=300)
plt.show()
