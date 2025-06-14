#!python3
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import chemgen as cg

# Load the reaction mechanism
gas = ct.Solution("FFCM2_Model.yaml")

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
C = 0
T = 0
for t in time:
    network.advance(t)
    temperature.append(reactor.T)
    data.append([t, reactor.T])
    if t>8e-6:
        T = reactor.T 
        C = reactor.thermo.concentrations
        break

big_dt = 2.5e-6
T_in = T 
C_in = C
Ts = []
Cs = []
Ts = []
Cs = []

for r in range(7):
    print(f"Refinement: {r+1}")
    dt = big_dt / 2**r
    n_steps = int(np.ceil(big_dt / dt))
    print(f"nsteps: {n_steps}")
    
    T = T_in
    C = C_in.copy()
    
    for k in range(n_steps):
        y = cg.sdirk4(C, T, dt, 1e-14, 10)
        C = y[1:]
        T = cg.temperature_from_internal_energy(C, y[0])
    
    Ts.append(T)
    Cs.append(C)

Ts = np.array(Ts)
Cs = np.array(Cs)
np.save("Tconvergence.npy",Ts)
np.save("Cconvergence.npy",Cs)
