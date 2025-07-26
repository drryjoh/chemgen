#!python3
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
try:
    plt.style.use('seaborn-colorblind')
except OSError:
    plt.style.use('seaborn-v0_8-colorblind')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
print(len(colors))

# Load the reaction mechanism
gas = ct.Solution("FFCM2_model.yaml")

# Define initial conditions
test_conditions = {
    "temperature": 1800,  # K
    "pressure": 101325.0,  # Pa
    "species": {
        "O2": 0.2,
        "N2": 0.4,
        "H2": 0.2,
        "C2H4": 0.2,
    }
}

gas.TPX = (
    test_conditions["temperature"],
    test_conditions["pressure"],
    test_conditions["species"]
)

# Create a reactor and insert the gas
reactor = ct.IdealGasMoleReactor(gas)
network = ct.ReactorNet([reactor])

# Define simulation time (in seconds)
time_end = 2e-5
dt_small = 1e-6
n_steps = int(time_end/dt_small)
#time_end = 200 * 2e-7  # Convert ns to seconds
#n_steps = 200  # Number of time steps
time = np.linspace(0, time_end, n_steps)

temperature = []
data = []
# Before the loop, set solver options
network.linear_solver_type = 'GMRES'           # Use GMRES instead of dense solver
network.preconditioner = ct.AdaptivePreconditioner()
network.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
network.rtol = 1e-10                           # Optional: tighter relative tolerance
network.atol = 1e-20                           # Optional: tighter absolute tolerance

# Now run the simulation loop

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

# Plot results
plt.plot(data[:, 0]*1000000.0, data[:, 1],'-', label = "Cantera", lw=5)
d = np.loadtxt("backward_euler.txt")
plt.plot(d[:, 0]*1000000.0, d[:, 1],'-d', label = "ChemGen Backward Euler", markevery=int(len(d[:, 0])/10))
d = np.loadtxt("sdirk2.txt")
plt.plot(d[:, 0]*1000000.0, d[:, 1],'-x', label = "ChemGen SDIRK2", markevery=int(len(d[:, 0])/10))
d = np.loadtxt("sdirk4.txt")
plt.plot(d[:, 0]*1000000.0, d[:, 1],':^', label = "ChemGen SDIRK4", markevery=int(len(d[:, 0])/10))
d = np.loadtxt("ros.txt")
plt.plot(d[:, 0]*1000000.0, d[:, 1],'--s',color='purple', label = "ChemGen Rosenbroc", markevery=int(len(d[:, 0])/10))

d = np.loadtxt("yass.txt")
plt.plot(d[:, 0]*1000000.0, d[:, 1],'-d',color='orange', label = "ChemGen YASS", markevery=int(len(d[:, 0])/10))

plt.legend()
plt.xlabel("Time ($\mu$s)")
plt.ylabel("Temperature (K)")
plt.title("Temperature Evolution in Homogeneous Reactor")
plt.savefig("implicit_time.png",dpi=300)
plt.show()
