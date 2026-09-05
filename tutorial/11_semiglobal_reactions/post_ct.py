#!python3
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Load the reaction mechanism
gas = ct.Solution("semiglobal_butadiene.yaml")

# Same initial condition as test_configuration.yaml / custom_test.py
gas.TPX = (
    1400,
    101325.0,
    {"C4H6": 0.1, "O2": 0.5, "CO": 0.1, "H2O": 0.1, "CO2": 0.2},
)

reactor = ct.IdealGasReactor(gas)
network = ct.ReactorNet([reactor])

# chem_out.txt holds both solvers' trajectories, tagged "RK4"/"SDIRK2" in column 0
rk4_rows = []
sdirk2_rows = []
with open("chem_out.txt") as f:
    for line in f:
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "RK4":
            rk4_rows.append([float(v) for v in parts[1:3]])
        elif parts[0] == "SDIRK2":
            sdirk2_rows.append([float(v) for v in parts[1:3]])
rk4 = np.array(rk4_rows)     # columns: t, T
sdirk2 = np.array(sdirk2_rows)

time_end = rk4[-1, 0]
time = np.linspace(0, time_end, 200)
temperature = []
for t in time:
    network.advance(t)
    temperature.append(reactor.T)

plt.plot(time * 1e6, temperature, '-r', label="Cantera")
plt.plot(rk4[:, 0] * 1e6, rk4[:, 1], '--k', label="ChemGen RK4")
plt.plot(sdirk2[:, 0] * 1e6, sdirk2[:, 1], ':bs', label="ChemGen SDIRK2", markevery=3)
plt.legend()
plt.xlabel("Time ($\\mu$s)")
plt.ylabel("Temperature (K)")
plt.title("Semi-global butadiene mechanism: homogeneous reactor")
plt.savefig("semiglobal_reactor.png", dpi=150)
plt.show()
