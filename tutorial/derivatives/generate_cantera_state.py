#!python3
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Load the reaction mechanism
gas = ct.Solution("ffcm2_h2.yaml")
gas = ct.Solution("burke.yaml")
gas = ct.Solution("h2o2.yaml")

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
n_steps = 300  # Number of time steps
time = np.linspace(0, time_end, n_steps)

temperature = []
data = []
printed = False
for t in time:
    network.advance(t)
    temperature.append(reactor.T)
    data.append([t, reactor.T] + list(reactor.thermo.concentrations))
    if reactor.T > 2400 and not printed:
        printed = True
        print(f"Pressure: {reactor.thermo.P}")
        print(f"Temperature: {reactor.T}")
        for k, specie in enumerate(gas.species_names):
            print(f"""    - name: {specie}
      MoleFraction: {reactor.thermo.X[k]}""")


    
data = np.array(data)
f = open("h2o2_data.csv","w")
for line in data:
    line_to_write = ','.join([str(s) for s in line])
    f.write(f"{line_to_write}\n")
f.close()
# Save results to file

# Plot results
#d = np.loadtxt("chem_out.txt")
plt.plot(data[:, 0]*1000.0, data[:, 1],'-r', label = "Cantera")
#plt.plot(d[:, 0]*1000.0, d[:, 1],'--k', label = "ChemGen")
#plt.legend()
plt.xlabel("Time ($\mu$s)")
plt.ylabel("Temperature (K)")
plt.title("Temperature Evolution in Homogeneous Reactor")
#plt.savefig("rk4.png")
#plt.show()
