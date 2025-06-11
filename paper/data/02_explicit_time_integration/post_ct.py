#!python3
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Load the reaction mechanism
gas = ct.Solution("ffcm2_h2.yaml")

# Define initial conditions
test_conditions = {
    "temperature": 1800,  # K
    "pressure": 101325.0,  # Pa
    "species": {
        "O2": 0.2,
        "N2": 0.6,
        "H2": 0.2
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
time_end = 200000 * 1e-9  # Convert ns to seconds
n_steps = 500  # Number of time steps
time = np.linspace(0, time_end, n_steps)

temperature = []
data = []

species_names = reactor.thermo.species_names  # optional: use if you want to label columns
for t in time:
    network.advance(t)
    temperature.append(reactor.T)
    concentrations = reactor.thermo.concentrations  # array of species concentrations [kmol/m^3]
    data.append([t, reactor.T] + list(concentrations))
data = np.array(data)

# Save results to file

# Plot results
d = np.loadtxt("chem_out.txt")
plt.plot(data[:, 0]*1000.0, data[:, 1],'-r', label = "Cantera")
plt.plot(d[:, 0]*1000.0, d[:, 1],'--k', label = "ChemGen")
plt.legend()
plt.xlim([0,0.2])
plt.xlabel("Time ($\mu$s)")
plt.ylabel("Temperature (K)")
plt.title("Temperature Evolution in Homogeneous Reactor")
plt.savefig("rk4.png")


#plot species
species_names = gas.species_names
species_index = {name: i for i, name in enumerate(species_names)}

fig, ax1 = plt.subplots()
time = data[:, 0]
concentration_data = data[:, 2:]  # skip time and temperature
color_idx = np.linspace(0,1,7)
for k, species in enumerate(["H2", "O2", "H2O"]):
    idx = species_index[species]
    ax1.plot(time*1000, concentration_data[:, idx], color = plt.cm.jet(color_idx[k]), label=species)
    ax1.plot(d[:, 0]*1000.0, d[:, 2+idx],'--k')

ax1.set_xlabel("Time ($\mu$s)")
ax1.set_ylabel('Major Species Concentration [kmol/m³]')
ax2 = ax1.twinx()
multipliers = [10000,10,10,10]
for k, (multiplier, species) in enumerate(zip(multipliers, ["H2O2", "OH", "H", "O"])):
    idx = species_index[species]
    ax2.plot(time*1000, multiplier * concentration_data[:, idx],  color = plt.cm.jet(color_idx[3 + k]), linestyle='-', label=f"{species} $\\times$ {multiplier}")
    if k == len(species) - 1:
        ax2.plot(d[:, 0]*1000.0, multiplier * d[:, 2+idx],'--k',label="ChemGen")
    else:
        ax2.plot(d[:, 0]*1000.0, multiplier * d[:, 2+idx],'--k')

ax1.set_xlim([0, 0.15])
ax2.set_xlim([0, 0.15])
ax2.set_ylabel('Minor Species Concentration [kmol/m³]')
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', ncol=3, fontsize=8)

plt.title('Species Concentrations Over Time')
plt.tight_layout()
plt.show()
