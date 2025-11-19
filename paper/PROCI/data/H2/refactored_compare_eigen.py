#!python3
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# Existing chemgen backends
import chemgen_conservative as cgc
import chemgen_temperature as cgt
import chemgen_dtemperature_ignore as cgc_it
import chemgen_conservative_ignore_species as cgc_is
import chemgen_temperature_ignore_species as cgt_is
# Example: when you add more, import and register here
# import chemgen_alt1 as cga
# import chemgen_alt2 as cgb

# ---------- Utilities ----------
def inv_largest_tau(J):
    # tau_min = 1 / max(|lambda|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(eig)
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size == 0 else 1.0 / np.max(nz)

def inv_smallest_tau(J):
    # tau_max = 1 / min(|lambda|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(eig)
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size == 0 else 1.0 / np.min(nz)

def stiffness(J):
    # kappa = tau_max / tau_min = max(|lambda|)/min(|lambda|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(eig)
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size < 2 else np.max(nz) / np.min(nz)

# Colorblind-friendly (Okabe–Ito)
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#000000", "#CC79A7", "#F0E442"]

# ---------- Register backends ----------
# Each backend supplies:
# - key: identifier
# - label: for legends
# - mod: imported module object
# - temp_from_state(y): extract temperature from state vector y returned by backwards_euler
# - style: plotting style

BACKENDS = [
    dict(
        key="cgc",
        label="Conservative",
        mod=cgc,
        temp_from_state=lambda y: cgc.temperature(y[1:], y[0]),
    ),
    dict(
        key="cgt",
        label="Temperature",
        mod=cgt,
        temp_from_state=lambda y: y[0],
    ),
    #dict(
    #    key="cgc_it",
    #    label="Conservative_nodT",
    #    mod=cgc_it,
    #    temp_from_state=lambda y: cgc_it.temperature(y[1:], y[0]),
    #),
    dict(
        key="cgc_is",
        label="Conservative_nodS",
        mod=cgc_is,
        temp_from_state=lambda y: cgc_it.temperature(y[1:], y[0]),
    ),
    dict(
        key="cgt_is",
        label="Temperature_nodS",
        mod=cgt_is,
        temp_from_state=lambda y: y[0],
    ),
    # When you add more, copy one of these blocks and adjust:
    # dict(key="cga", label="Alt-1", mod=cga, temp_from_state=lambda y: ...),
    # dict(key="cgb", label="Alt-2", mod=cgb, temp_from_state=lambda y: ...),
]

# Assign styles deterministically
for i, be in enumerate(BACKENDS):
    be["style"] = dict(color=PALETTE[i % len(PALETTE)],
                       ls="--" if i % 2 else "-.",
                       lw = 2)

# Choose a baseline for ratio plots (defaults to first backend if not found)
BASELINE_KEY = "cgc" if any(b["key"] == "cgc" for b in BACKENDS) else BACKENDS[0]["key"]

# ---------- Cantera setup ----------
gas = ct.Solution("ffcm2_h2.yaml")
gas.set_equivalence_ratio(phi=1.0, fuel="H2", oxidizer="O2:1.0, N2:3.76")  # O2:N2 = 1:3.76
gas.TP = 2500, 101325
reactor = ct.IdealGasReactor(gas)
net = ct.ReactorNet([reactor])

n_steps = 800
dt = 1e-7
time = np.linspace(0.0, n_steps * dt, n_steps)

# Reference trajectory from Cantera
T_phys = []

# Per-backend storage
store = {}
for be in BACKENDS:
    store[be["key"]] = dict(
        T_pred=[],
        inv_tau_min=[],
        inv_tau_max=[],
        kappa=[],
        lin_iter=[],
        newt_iter=[],
    )

# ---------- Time integration loop ----------
for t in time:
    T = reactor.T
    C = reactor.thermo.concentrations

    # Run each backend at the current physical state
    for be in BACKENDS:
        mod = be["mod"]
        y, n_lin, n_newt = mod.backwards_euler(C, T, dt)
        T_hat = be["temp_from_state"](y)

        s = store[be["key"]]
        s["T_pred"].append(T_hat)
        s["lin_iter"].append(n_lin)
        s["newt_iter"].append(n_newt)

    # Advance the physical reactor to this time
    net.advance(t)
    T_phys.append(reactor.T)

    # Form Jacobians at the updated physical state
    T_now = reactor.T
    C_now = reactor.thermo.concentrations
    for be in BACKENDS:
        J = np.array(be["mod"].source_jacobian(C_now, T_now))
        s = store[be["key"]]
        s["inv_tau_min"].append(inv_largest_tau(J))
        s["inv_tau_max"].append(inv_smallest_tau(J))
        s["kappa"].append(stiffness(J))

# Convert to arrays
T_phys = np.array(T_phys)
for be in BACKENDS:
    s = store[be["key"]]
    for k in s:
        s[k] = np.array(s[k])

# ---------- Plots ----------
fig, axs = plt.subplots(4, 1, figsize=(10, 7), sharex=True)

# 1) Temperature
axs[0].plot(time, T_phys, color="#000000", label="Cantera")
for be in BACKENDS:
    st = be["style"]
    axs[0].plot(time, store[be["key"]]["T_pred"], label=f"{be['label']}", **st)
axs[0].set_ylabel("Temperature [K]")
axs[0].set_title("Temperature evolution")
axs[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 2 )

# 2) Inverse eigenvalues: tau_min and tau_max per backend
for be in BACKENDS:
    st = be["style"]
    axs[1].semilogy(time, store[be["key"]]["inv_tau_min"],
                    label=rf"$\tau_{{min}}$ {be['label']}", **st)
    # derive variant style for tau_max: dotted of the same color
    axs[1].semilogy(time, store[be["key"]]["inv_tau_max"],
                    linestyle=":", color=st["color"],
                    label=rf"$\tau_{{max}}$ {be['label']}")
axs[1].set_ylabel(r"Characteristic time $\tau = 1/\lambda$")
axs[1].set_title("Inverse of largest and smallest nonzero eigenvalues")
axs[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 2)

# 3) Stiffness per backend
for be in BACKENDS:
    st = be["style"]
    axs[2].semilogy(time, store[be["key"]]["kappa"],
                    label=f"{be['label']}", **st)
axs[2].set_ylabel(r"Stiffness, $\tau_{max}/\tau_{min}$")
axs[2].set_xlabel("Time [s]")
axs[2].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 2)

# 4) Stiffness ratios vs baseline
baseline = store[BASELINE_KEY]["kappa"]
for be in BACKENDS:
    st = be["style"]
    ratio = store[be["key"]]["kappa"] / baseline
    axs[3].semilogy(time, ratio, label=f"{be['label']} / {BASELINE_KEY}", **st)
axs[3].set_ylabel("Stiffness ratio\nvs baseline")
axs[3].set_xlabel("Time [s]")
axs[3].set_ylim([0.1, np.nanmax([np.nanmax(store[b['key']]['kappa']/baseline) for b in BACKENDS]) * 10.0])

# Grids and reference line at 10^0 = 1
for ax in axs:
    ax.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
    ax.grid(True, which="minor", linestyle="--", linewidth=0.3, alpha=0.5)
axs[3].axhline(1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.8)

plt.tight_layout()
plt.subplots_adjust(right=0.78)
plt.savefig("stiffness_formulation.png", dpi=300, bbox_inches="tight")

# Iteration plots
fig, axs = plt.subplots(3, 1, figsize=(10, 7), sharex=True)

# Linear iterations per step
for be in BACKENDS:
    st = be["style"]
    total = np.nansum(store[be["key"]]["lin_iter"])
    axs[0].plot(time, store[be["key"]]["lin_iter"],
                label=f"{be['label']} total = {int(total)}", **st)
axs[0].set_ylabel("Linear iterations\nper time step")
axs[0].set_xlabel("Time [s]")
axs[0].legend(ncol=2)

# Nonlinear iterations per step
for be in BACKENDS:
    st = be["style"]
    axs[1].plot(time, store[be["key"]]["newt_iter"],
                label=f"{be['label']}", **st)
axs[1].set_ylabel("Nonlinear iterations\nper time step")
axs[1].set_xlabel("Time [s]")
axs[1].legend(ncol=2)

# Linear per nonlinear
for be in BACKENDS:
    st = be["style"]
    num = store[be["key"]]["lin_iter"].astype(float)
    den = store[be["key"]]["newt_iter"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(den > 0, num / den, np.nan)
    axs[2].plot(time, ratio, label=f"{be['label']}", **st)
axs[2].set_ylabel("Linear iterations\nper nonlinear step")
axs[2].set_xlabel("Time [s]")
axs[2].legend(ncol=2)

plt.tight_layout()
plt.show()
