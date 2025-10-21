#!python3
import chemgen_conservative as cgc
import chemgen_temperature as cgt
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
def process_large(J):
    eigvals = np.linalg.eigvals(J)
    mag = np.abs(eigvals)
    tol = 1e-12 * np.max(mag)
    nonzero = mag[mag > tol]
    return 1/np.max(nonzero)
def process_small(J):
    eigvals = np.linalg.eigvals(J)
    mag = np.abs(eigvals)
    tol = 1e-12 * np.max(mag)
    nonzero = mag[mag > tol]
    return 1/np.min(nonzero)
def stiffness(J):
    eigvals = np.linalg.eigvals(J)
    mag = np.abs(eigvals)
    tol = 1e-12 * np.max(mag)
    nonzero = mag[mag > tol]
    stiffness_ratio = np.max(nonzero) / np.min(nonzero)
    return stiffness_ratio

import cantera as ct
gas = ct.Solution("ffcm2_h2.yaml")
gas.set_equivalence_ratio(phi=1.0, fuel="H2", oxidizer="O2:1, N2:3.76")
gas.TP = 2200, 101325  # temperature [K], pressure [Pa]
reactor = ct.IdealGasReactor(gas)

network = ct.ReactorNet([reactor])

n_steps = 400
dt = 2e-7
time_end = n_steps * dt

time = np.linspace(0, time_end, n_steps)

temperature = []
temperature_cgc = []
temperature_cgt = []
inv_largest_eigs_c = []
inv_largest_eigs_t = []
inv_smallest_eigs_c = []
inv_smallest_eigs_t = []
cond_c = []
cond_t = []
iterations_cgc = []
iterations_cgt = []
newt_iterations_cgc = []
newt_iterations_cgt = []

for t in time:
    T = reactor.T
    C = reactor.thermo.concentrations
    y_cgc, n_cgc, n_cgc_newt  = cgc.backwards_euler(C, T, dt)
    y_cnt, n_cgt, n_cgt_newt = cgt.backwards_euler(C, T, dt)
    iterations_cgc.append(n_cgc)
    iterations_cgt.append(n_cgt)
    newt_iterations_cgc.append(n_cgc_newt)
    newt_iterations_cgt.append(n_cgt_newt)

    temperature_cgt.append(y_cnt[0])
    temperature_cgc.append(cgc.temperature(y_cgc[1:],y_cgc[0]))
    network.advance(t)
    T = reactor.T
    temperature.append(T)

    Jc = np.array(cgc.source_jacobian(C, T))
    Jt = np.array(cgt.source_jacobian(C, T))
    inv_largest_eigs_c.append(process_large(Jc))
    inv_largest_eigs_t.append(process_large(Jt))
    inv_smallest_eigs_c.append(process_small(Jc))
    inv_smallest_eigs_t.append(process_small(Jt))
    cond_c.append(stiffness(Jc))
    cond_t.append(stiffness(Jt))


temperature = np.array(temperature)
temperature_cgt = np.array(temperature_cgt)
temperature_cgc = np.array(temperature_cgc)

iterations_cgc = np.array(iterations_cgc)
iterations_cgt = np.array(iterations_cgt)
newt_iterations_cgt = np.array(newt_iterations_cgt)
newt_iterations_cgc = np.array(newt_iterations_cgc)

inv_largest_eigs_c = np.array(inv_largest_eigs_c)
inv_largest_eigs_t = np.array(inv_largest_eigs_t)
inv_smallest_eigs_c = np.array(inv_smallest_eigs_c)
inv_smallest_eigs_t = np.array(inv_smallest_eigs_t)
cond_c  = np.array(cond_c)
cond_t  = np.array(cond_t)

# Colorblind-friendly palette (Okabe–Ito)
colors = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "sky": "#56B4E9",
    "black": "#000000"
}

# Plot results
fig, axs = plt.subplots(4, 1, figsize=(10, 7), sharex=True)

# --- Temperature ---
axs[0].plot(time, temperature, color=colors["black"])
axs[0].plot(time, temperature_cgt, color=colors["orange"],linestyle = '--')
axs[0].plot(time, temperature_cgc, color=colors["red"], linestyle = '-.')
axs[0].set_ylabel("Temperature [K]")
axs[0].set_title("Temperature evolution")

# --- Inverse eigenvalues ---
axs[1].semilogy(time, inv_largest_eigs_c, color=colors["blue"], linestyle='-', label=r"$\tau_{min}$ cons-form")
axs[1].semilogy(time, inv_largest_eigs_t, color=colors["orange"], linestyle='--', label=r"$\tau_{min}$ T-form")
axs[1].semilogy(time, inv_smallest_eigs_c, color=colors["green"], linestyle='-', label=r"$\tau_{max}$ cons-form")
axs[1].semilogy(time, inv_smallest_eigs_t, color=colors["red"], linestyle='--', label=r"$\tau_{max}$ T-form")
axs[1].set_ylabel(r"Characteristic time $\tau = 1/\lambda$")
axs[1].set_title("Inverse of largest and smallest non-zero eigenvalues")

# --- Condition number / stiffness ---
axs[2].semilogy(time, cond_c, color=colors["blue"], linestyle='-', label="Cons-form")
axs[2].semilogy(time, cond_t, color=colors["orange"], linestyle='--', label="T-form")
axs[2].set_ylabel("Stiffness,\n$\\tau_{max}/\\tau_{min}$")
axs[2].set_xlabel("Time [s]")

# --- Condition number / stiffness ratio ---
axs[3].semilogy(time, cond_t / cond_c,
                color=colors["blue"],
                linestyle='-',
                label="Ratio of stiffness")

axs[3].set_ylabel(r"Ratio of Stiffness,"
                  "\n"
                  r"$\frac{\tau_{\max}}{\tau_{\min}}_{T} / \frac{\tau_{\max}}{\tau_{\min}}_{C}$")
axs[3].set_xlabel("Time [s]")
axs[3].set_ylim([0.1, np.max(cond_t / cond_c) * 10])

# Enable both major and minor gridlines for better readability
axs[3].grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.7)
axs[3].grid(True, which='minor', linestyle='--', linewidth=0.3, alpha=0.5)
axs[3].axhline(1.0, color='gray', linestyle=':', linewidth=1.0, alpha=0.8)


# Move legends outside (right side)
axs[1].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
axs[2].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

# Layout and save
plt.tight_layout()
plt.subplots_adjust(right=0.8)
plt.savefig("stiffness_formulation.png", dpi=300, bbox_inches='tight')

fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
axs[0].plot(time, iterations_cgt, color=colors["orange"],linestyle = '--', label = f"Temperature Formulation Total = {np.sum(iterations_cgt)}")
axs[0].plot(time, iterations_cgc, color=colors["red"], linestyle = '-.', label = f"Conservative Formulation Total =  {np.sum(iterations_cgc)}")
axs[0].set_ylabel("Linear iterations\n per time step")
axs[0].set_xlabel("Time [s]")
axs[0].legend()

axs[1].plot(time, newt_iterations_cgt, color=colors["orange"],linestyle = '--', label = "Temperature Formulation")
axs[1].plot(time, newt_iterations_cgc, color=colors["red"], linestyle = '-.', label = "Conservative Formulation")
axs[1].set_ylabel("Nonlinear iterations\n per time step")
axs[1].set_xlabel("Time [s]")
axs[1].legend()

axs[2].plot(time, iterations_cgt/newt_iterations_cgt, color=colors["orange"],linestyle = '--', label = "Temperature Formulation")
axs[2].plot(time, iterations_cgc/newt_iterations_cgc, color=colors["red"], linestyle = '-.', label = "Conservative Formulation")
axs[2].set_ylabel("Linear iterations\n per nonlinear step")
axs[2].set_xlabel("Time [s]")
axs[2].legend()
plt.show()
