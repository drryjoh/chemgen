#!python3
import chemgen as cg
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution("ffcm2_h2.yaml")
gas.TPX  = 1800, 1013250, "H2:0.2, O2:0.1, N2:0.7"
reactor = ct.IdealGasReactor(gas)
network = ct.ReactorNet([reactor])

n_steps = 400
time_end = n_steps * 2e-8

time = np.linspace(0, time_end, n_steps)

temperature = []
inv_largest_eigs = []

for t in time:
    network.advance(t)
    T = reactor.T
    C = reactor.thermo.concentrations
    temperature.append(T)

    J = np.array(cg.source_jacobian(C, T))
    eigs = np.linalg.eigvals(J)

    # Filter out near-zero eigenvalues
    non_zero_eigs = [e.real for e in eigs if np.abs(e) > 1e-12]
    if non_zero_eigs:
        largest = max(np.abs(non_zero_eigs))
        smallest = min(np.abs(non_zero_eigs))
        inv_largest_eigs.append(1/largest)
    else:
        inv_largest_eigs.append(np.nan)  # or 0.0

temperature = np.array(temperature)
inv_largest_eigs = np.array(inv_largest_eigs)

# Plot results
plt.figure(figsize=(10, 5))

plt.subplot(2, 1, 1)
plt.plot(time, temperature, '-k')
plt.ylabel("Temperature [K]")
plt.title("Temperature evolution")

plt.subplot(2, 1, 2)
plt.semilogy(time, inv_largest_eigs, '-r')
plt.ylabel("1 / λ_max")
plt.xlabel("Time [s]")
plt.title("Inverse of largest non-zero eigenvalue")

plt.tight_layout()
plt.show()
