#!python3
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Load the reaction mechanism
gas = ct.Solution("ffcm2_h2.yaml")
plt.rcParams["lines.linewidth"] = 2.5


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
n_steps = 500
time = np.linspace(0, time_end, n_steps)

temperature = []
data = []

species_names = reactor.thermo.species_names
for t in time:
    network.advance(t)
    temperature.append(reactor.T)
    concentrations = reactor.thermo.concentrations
    data.append([t, reactor.T] + list(concentrations))
data = np.array(data)

# Load ChemGen output
d = np.loadtxt("chem_out.txt")
species_index = {name: i for i, name in enumerate(species_names)}
concentration_data = data[:, 2:]  # exclude time and temperature

# Set up side-by-side plots
fig, (ax_temp, ax_species) = plt.subplots(1, 2, figsize=(12, 5))

# ---- Plot Temperature Evolution (Left) ----
ax_temp.plot(data[:, 0]*1000.0, data[:, 1], '-r', label="Cantera")
ax_temp.plot(d[:, 0]*1000.0, d[:, 1], '--k', label="ChemGen")
np.save("c_T.npy", data[:, 1])
np.save("c_time.npy", data[:, 0])
ax_temp.set_xlim([0, 0.2])
ax_temp.set_xlabel("Time ($\mu$s)")
ax_temp.set_ylabel("Temperature (K)")
ax_temp.set_title("Temperature Evolution")
ax_temp.legend()

# ---- Plot Species Concentrations (Right) ----
color_idx = np.linspace(0, 1, 7)
# Major species on left y-axis
for k, species in enumerate(["H2", "O2", "H2O"]):
    idx = species_index[species]
    ax_species.plot(time*1000, concentration_data[:, idx], 
                    color=plt.cm.tab10(color_idx[k]), label=species)
    ax_species.plot(d[:, 0]*1000.0, d[:,2 + idx], '--k')

ax_species.set_xlabel("Time ($\mu$s)")
ax_species.set_ylabel('Major Species Concentration [kmol/m³]')
ax_species.set_xlim([0, 0.125])

# Minor species on right y-axis
ax2 = ax_species.twinx()
multipliers = [10000, 10, 10, 10]
for k, (multiplier, species) in enumerate(zip(multipliers, ["H2O2", "OH", "H", "O"])):
    idx = species_index[species]
    ax2.plot(time*1000, multiplier * concentration_data[:, idx],
             color=plt.cm.tab10(color_idx[3 + k]), linestyle='-', label=f"{species} $\\times$ {multiplier}")
    ax2.plot(d[:, 0]*1000.0, multiplier * d[:, 2+idx], '--k')

ax2.set_ylabel('Minor Species Concentration [kmol/m³]')
ax2.set_xlim([0, 0.125])

# Combine legends
lines1, labels1 = ax_species.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax_species.legend(lines1 + lines2, labels1 + labels2, loc='upper right', ncol=2, fontsize=8)

ax_species.set_title('Species Concentrations')

# Final layout adjustments
plt.tight_layout()
plt.savefig("rk4.png", dpi=300)
plt.show()
