import cantera as ct

# Create a solution object from the given chem.yaml file
solution = ct.Solution("chemkin/chem.yaml")

# Define the initial conditions
temperature = 1000  # K
pressure = 202650  # Pa
mole_fractions = "H2:1, N2:3.76, O2:1"

# Set the state of the solution
solution.TPY = 2477.41, 460544, [0.000113656,1.55054e-05,0.11162,0.00109796,0.00661249,2.00938e-05,1.06531e-06,0.124623,0,0.755896]

#temperature, pressure, mole_fractions

# Create a constant-volume reactor
reactor = ct.Reactor(solution)

# Create a reactor network
network = ct.ReactorNet([reactor])

# Define simulation time parameters
time = 0.0  # start time [s]
end_time = 0.000333333  # end time [s]
deltaT = 0.000333333
data = []  # list to store time and state information

# Integrate in time
c0 = reactor.thermo.concentrations
print("concentrations")
print(c0)
print(reactor.thermo.density*reactor.thermo.Y/reactor.thermo.molecular_weights)
while time < end_time:
    
    network.advance(time + deltaT)  # advance the simulation by a fixed time step
    time = time + deltaT
    data.append((time, reactor.T, reactor.thermo.P, reactor.thermo['H2'].X[0], (reactor.thermo.concentrations-c0)/(deltaT)*reactor.thermo.molecular_weights))  # store time, temperature, pressure, and H2 mole fraction

# Print final results
print("Final State:")
print(f"Time: {time:.6f} s")
print(f"Temperature: {reactor.T:.2f} K")
print(f"Pressure: {reactor.thermo.P:.2f} Pa")
print(
"""
H2: -0.00127918, 0.000334075, -0.261367
H: 0.0150064, 6.49144e-05, 0.00432548
O2: 0.267808, -0.00594328, -0.0221922
O: -0.24088, 0.00330528, -0.0137217
OH: 0.461698, 0.0122945, 0.0266287
HO2: -0.125944, 2.89354e-05, -0.000229749
H2O2: -0.092654, 1.3139e-06, -1.41808e-05
H2O: -0.283755, -0.0100857, 0.0355438
""")
# Optionally, save the results to a file
with open("reactor_results.csv", "w") as file:
    file.write("time,temperature,pressure,H2_mole_fraction\n")
    for t, T, P, X_H2, RR in data:
        file.write(f"{t:.6e},{T:.2f},{P:.2f},{X_H2:.6e}\n")
        print(RR)

