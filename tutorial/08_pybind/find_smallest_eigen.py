#!python3
import chemgen as cg
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

def find_eigs(J):
    # Filter out near-zero eigenvalues
    eigs = np.linalg.eigvals(J)
    non_zero_eigs = [e.real for e in eigs if np.abs(e) > 1e-12]
    if non_zero_eigs:
        largest = max(np.abs(non_zero_eigs))
        smallest = min(np.abs(non_zero_eigs))
        return 1/largest
    else:
        return np.nan  # or 0.0
    return 0.0

gas = ct.Solution("andrea.yaml")
gas.TPX  = 1200, 101325, "CH4:0.01, O2:0.02, AR:0.97"
reactor = ct.IdealGasReactor(gas)
network = ct.ReactorNet([reactor])

time_end = 1.5
n_steps = int(1000)
dt = time_end/n_steps
time_end = n_steps * dt

time = np.linspace(0, time_end, n_steps)

temperature = []
fastest_reactions = []
fastest_times = []
fastest_times_removed = []
for t in time:
    network.advance(t)
    T = reactor.T
    C = reactor.thermo.concentrations
    temperature.append(T)

    Jtotal = np.array(cg.source_jacobian(C, T))
    smallest_time_scale_total = find_eigs(Jtotal)

    #fastest_reaction_at_time = [smallest_time_scale_total]  # include total baseline
    fastest_reaction_at_time = []
    fastest_times.append(smallest_time_scale_total)
    for i in range(gas.n_reactions):
        J = np.array(cg.source_jacobian_remove_reaction(C, T, Jtotal, i))
        time_at_i = find_eigs(J)
        if (time_at_i - smallest_time_scale_total)/smallest_time_scale_total > 0.1:
            fastest_reaction_at_time.append(i)

    fastest_times_removed.append(fastest_reaction_at_time)

temperature = np.array(temperature)
fastest_times = np.array(fastest_times)

# Plot results
plt.figure(figsize=(10, 5))

plt.subplot(2, 1, 1)
plt.plot(time, temperature, '-k')
plt.ylabel("Temperature [K]")
plt.title("Temperature evolution")

plt.subplot(2, 1, 2)
for i, fastest_reactions_at_time in enumerate(fastest_times_removed):
    for reaction in fastest_reactions_at_time:
        plt.plot(time[i], reaction, 'ok')
#plt.semilogy(time, fastest_times,'-r')
plt.ylabel("Fastest Reactions")
plt.xlabel("Time [s]")
plt.title("Inverse of largest non-zero eigenvalue")

plt.tight_layout()
plt.show()
