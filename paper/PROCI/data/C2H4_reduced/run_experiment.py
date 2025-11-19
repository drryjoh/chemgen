#!python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import cantera as ct

# User chemgen backends
import chemgen_conservative as cgc
import chemgen_temperature as cgt
# Optional imports (kept commented)
# import chemgen_dtemperature_ignore as cgc_it
# import chemgen_conservative_ignore_species as cgc_is
# import chemgen_temperature_ignore_species as cgt_is

from numerics import backwards_euler, sdirk2, jacobi_precondition, full_jacobian
from metrics import (inv_largest_tau, inv_smallest_tau, eigs, svs, alphas, stiffness,
                     departure, clustering_scores, PALETTE)
from matplotlib.patches import Ellipse

import numpy as np

def draw_density_ellipse_curve(ax, xdata, ydata, color='red', lw=1.5, scale=2.0, npts=256, **plot_kwargs):
    """
    Draws a covariance ellipse (≈ density region) for log-log data by plotting a curve.
    The covariance is computed in log10 space and mapped back to data space.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    xdata, ydata : array-like
        Data arrays; only positive values are used.
    color : str
        Line color.
    lw : float
        Line width.
    scale : float
        # of std devs ~ 2 for ~95%, 3 for ~99%.
    npts : int
        Number of points along the ellipse.
    **plot_kwargs : passed to ax.plot
    """
    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)

    mask = (xdata > 0) & (ydata > 0) & np.isfinite(xdata) & np.isfinite(ydata)
    if mask.sum() < 3:
        return None  # not enough points

    X = np.log10(xdata[mask])
    Y = np.log10(ydata[mask])

    mean = np.array([X.mean(), Y.mean()])
    cov = np.cov(X, Y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]

    # Radii in log space
    rx, ry = scale * np.sqrt(np.maximum(vals, 0.0))

    # Parametric ellipse in log space
    t = np.linspace(0, 2*np.pi, npts, endpoint=True)
    R = vecs @ np.diag([rx, ry])
    pts_log = mean[:, None] + R @ np.vstack((np.cos(t), np.sin(t)))

    # Map back to data space
    xl = 10**pts_log[0, :]
    yl = 10**pts_log[1, :]

    line, = ax.plot(xl, yl, color=color, lw=lw, **plot_kwargs)
    return line



def clustering_scores_log(x, alpha=0.3):
    x = np.asarray(x, dtype=float)
    eps = 1e-30
    x_log = np.log10(np.abs(x) + eps)
    x_log = x_log[np.isfinite(x_log)]
    if len(x_log) < 2:
        return {"R_NN": np.nan, "CV_NN": np.nan, "F_eps": np.nan}
    x_log = np.sort(x_log)
    d = np.empty_like(x_log)
    d[0]  = x_log[1] - x_log[0]
    d[-1] = x_log[-1] - x_log[-2]
    d[1:-1] = np.minimum(x_log[1:-1]-x_log[:-2], x_log[2:]-x_log[1:-1])
    mean_d = np.mean(d)
    cv_nn  = np.std(d) / mean_d
    norm   = (x_log[-1]-x_log[0]) / (len(x_log)-1)
    R_nn   = mean_d / norm
    F_eps  = np.mean(d < alpha * norm)
    return {"R_NN": R_nn, "CV_NN": cv_nn, "F_eps": F_eps}


# ---------- Register backends ----------
BACKENDS = [
    dict(
        key="cgc",
        label="Conservative",
        require_internal_energy="yes",
        mod=cgc,
        temp_from_state=lambda y: cgc.temperature(y[1:], y[0]),
    ),
    dict(
        key="cgt",
        label="Temperature",
        require_internal_energy="no",
        mod=cgt,
        temp_from_state=lambda y: y[0],
    ),
]

# Assign styles deterministically
for i, be in enumerate(BACKENDS):
    be["style"] = dict(color=PALETTE[i % len(PALETTE)],
                       ls="--" if i % 2 else "-.",
                       lw=2)

BASELINE_KEY = "cgc" if any(b["key"] == "cgc" for b in BACKENDS) else BACKENDS[0]["key"]

# ---------- Cantera setup ----------
gas = ct.Solution("mechanism.yaml")
gas.set_equivalence_ratio(phi=1.0, fuel="C2H4", oxidizer="O2:1.0, N2:3.76")
gas.TP = 2500, 101325
reactor = ct.IdealGasReactor(gas)
net = ct.ReactorNet([reactor])

n_steps = 200
dt = 2e-7
time = np.linspace(0.0, n_steps * dt, n_steps)

# Reference trajectory from Cantera
T_phys = []

# Per-backend storage
store = {be["key"]: dict(T_pred=[], inv_tau_min=[], inv_tau_max=[], kappa=[],
                         lin_iter=[], newt_iter=[], eigs=[], svs=[], departure=[], alphas=[])
         for be in BACKENDS}

# ---------- Time integration loop ----------
for t in time:
    T = reactor.T
    C = reactor.thermo.concentrations

    for be in BACKENDS:
        # Choose integrator
        # y, n_lin, n_newt = sdirk2(C, T, dt, be, gmres_method="numpy")
        y, n_lin, n_newt = backwards_euler(C, T, dt, be, gmres_method="eigen")
        T_hat = be["temp_from_state"](y)

        s = store[be["key"]]
        s["T_pred"].append(T_hat)
        s["lin_iter"].append(n_lin)
        s["newt_iter"].append(n_newt)

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
        s["eigs"].append(eigs(full_jacobian(J, dt)))
        s["svs"].append(svs(full_jacobian(J, dt)))
        s["alphas"].append(alphas(full_jacobian(J, dt), be["mod"].source(C_now, T_now)))
        s["kappa"].append(stiffness(J))
        s["departure"].append(departure(J))

# Convert to arrays
T_phys = np.array(T_phys)
for be in BACKENDS:
    s = store[be["key"]]
    for k in s:
        s[k] = np.array(s[k])

# For plotting: microseconds
time_us = time * 1e6


# ---------- Singular values vs time + clustering metrics ----------
fig, ax = plt.subplots(2, 2, figsize=(15, 10), sharey=False)
#axs = [axs]
axs = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]
#axs = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]

for be in BACKENDS:
    if be["label"] in ("Conservative", "Temperature"):
        idx = 0 if be["label"] == "Conservative" else 2
        #marker = "+" if be["label"] == "Conservative" else "o"
        marker = "o" if be["label"] == "Conservative" else "o"
        key = be["key"]
        svs_all = store[key]["svs"]
        st = be["style"]
        for i, t_us in enumerate(time_us):
            svs_t = np.array(svs_all[i])
            # Windowing as in user code
            for _ in range(12):
                mean_log = np.mean(np.log(svs_t))
                idx_rem = np.argmax(np.abs(np.log(svs_t) - mean_log))
                svs_t = np.delete(svs_t, idx_rem)

            if svs_t.size == 0:
                continue
            axs[idx].semilogy(svs_t*0 + t_us, svs_t, marker,
                              mfc=mcolors.to_rgba(st["color"], alpha=0.7),
                              mec='black', markersize=4, label=be["label"] if i == 0 else None)
            cluster = clustering_scores_log(svs_t)
            axs[1].plot(t_us, cluster["R_NN"], 'o',
                        mfc=mcolors.to_rgba(st["color"], alpha=0.7),
                        mec='black', markersize=4, label=be["label"] if i == 0 else None)
            axs[1].set_ylabel("std(singular values)")
            axs[3].plot(t_us, cluster["CV_NN"], 'o',
                        mfc=mcolors.to_rgba(st["color"], alpha=0.7),
                        mec='black', markersize=4, label=be["label"] if i == 0 else None)
            axs[3].set_ylabel("CV_NN")

for a in axs:
    a.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
    a.grid(True, which="minor", linestyle="--", linewidth=0.3, alpha=0.5)

plt.tight_layout()
plt.savefig("svs_clustering.png", dpi=200, bbox_inches="tight")

# ---------- Singular values vs time + clustering metrics ----------
fig, ax = plt.subplots(1, 1, figsize=(15, 10), sharey=False)
axs = [ax]
#axs = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]

for be in BACKENDS:
    if be["label"] in ("Conservative", "Temperature"):
        idx = 0 if be["label"] == "Conservative" else 2
        #marker = "+" if be["label"] == "Conservative" else "o"
        marker = "o" if be["label"] == "Conservative" else "o"
        key = be["key"]
        svs_all = store[key]["svs"]
        eigs_all = store[key]["eigs"]
        alphas_all = store[key]["alphas"]
        np.save(f"{be['label']}_svs.npy",eigs_all)
        np.save(f"{be['label']}_alphas.npy",alphas_all)
        #mask = (svs_all >= 1e2) & (svs_all <= 1e12)
        #svs_all, alphas_all = svs_all[mask], alphas_all[mask]
        st = be["style"]
        axs[0].loglog(eigs_all, alphas_all, marker, mfc=mcolors.to_rgba(st["color"], alpha=0.7), mec='black', markersize=4)
        axs[0].set_xlabel("Singular Values, $\sigma_i$")
        axs[0].set_ylabel("$b$-alignment, $\\alpha_i=U^Tb$")
        draw_density_ellipse_curve(axs[0], svs_all, alphas_all, color=st["color"], scale=4.0, lw=4.0)



# ---------- Temperature and stiffness plots ----------
fig, axs = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

# 1) Temperature
axs[0].plot(time_us, T_phys, color="#000000", label="Cantera")
for be in BACKENDS:
    st = be["style"]
    axs[0].plot(time_us, store[be["key"]]["T_pred"], label=f"{be['label']}", **st)
axs[0].set_ylabel("Temperature [K]")
axs[0].set_title("Temperature evolution")
axs[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=2, fontsize=10)

# 2) Inverse eigenvalues: tau_min and tau_max per backend
for be in BACKENDS:
    st = be["style"]
    axs[1].semilogy(time_us, store[be["key"]]["inv_tau_min"],
                    label=rf"$\tau_{{min}}$ {be['label']}", **st)
    axs[1].semilogy(time_us, store[be["key"]]["inv_tau_max"],
                    linestyle=":", color=st["color"],
                    label=rf"$\tau_{{max}}$ {be['label']}")
axs[1].set_ylabel(r"Characteristic time $\tau = 1/\lambda$")
axs[1].set_title("Inverse of largest and smallest nonzero eigenvalues")
axs[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=2, fontsize=10)

# 3) Stiffness per backend
for be in BACKENDS:
    st = be["style"]
    axs[2].semilogy(time_us, store[be["key"]]["kappa"],
                    label=f"{be['label']}", **st)
axs[2].set_ylabel(r"Stiffness, $\tau_{max}/\tau_{min}$")
axs[2].set_xlabel(r"Time [$\mu$s]")
axs[2].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=2)

# 4) Nonnormality proxy
for be in BACKENDS:
    st = be["style"]
    ratio = store[be["key"]]["departure"]
    axs[3].plot(time_us, ratio, label=f"{be['label']}", **st)
axs[3].set_ylabel("Nonnormality proxy")
axs[3].set_xlabel(r"Time [$\mu$s]")
axs[3].legend()

for ax in axs:
    ax.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
    ax.grid(True, which="minor", linestyle="--", linewidth=0.3, alpha=0.5)
axs[3].axhline(1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.8)

plt.tight_layout()
plt.subplots_adjust(right=0.78)
plt.savefig("stiffness_formulation.png", dpi=300, bbox_inches="tight")

# ---------- Iteration plots ----------
fig, axs = plt.subplots(3, 1, figsize=(15, 10), sharex=True)

for be in BACKENDS:
    st = be["style"]
    np.save(f"time_{be['key']}.npy", time_us)
    np.save(f"linear_iterations_{be['key']}.npy", store[be["key"]]["lin_iter"])
    total = np.nansum(store[be["key"]]["lin_iter"])
    axs[0].plot(time_us, store[be["key"]]["lin_iter"],
                label=f"{be['label']} total = {int(total)}", **st)
axs[0].set_ylabel("Linear iterations\nper time step")
axs[0].set_xlabel(r"Time [$\mu$s]")
axs[0].legend(ncol=2)

for be in BACKENDS:
    st = be["style"]
    axs[1].plot(time_us, store[be["key"]]["newt_iter"],
                label=f"{be['label']}", **st)
axs[1].set_ylabel("Nonlinear iterations\nper time step")
axs[1].set_xlabel(r"Time [$\mu$s]")
axs[1].legend(ncol=2)

for be in BACKENDS:
    st = be["style"]
    num = store[be["key"]]["lin_iter"].astype(float)
    den = store[be["key"]]["newt_iter"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(den > 0, num / den, np.nan)
    axs[2].plot(time_us, ratio, label=f"{be['label']}", **st)
axs[2].set_ylabel("Linear iterations\nper nonlinear step")
axs[2].set_xlabel(r"Time [$\mu$s]")
axs[2].legend(ncol=2)

plt.tight_layout()
plt.savefig("iteration_metrics.png", dpi=300, bbox_inches="tight")
print("Finished. Artifacts: eigs_hist.png, svs_clustering.png, stiffness_formulation.png, iteration_metrics.png")

# ---------- Plots ----------
fig, axs = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

# 1) Temperature
axs[0].plot(time, T_phys, color="#000000", label="Cantera")
for be in BACKENDS:
    st = be["style"]
    axs[0].plot(time, store[be["key"]]["T_pred"], label=f"{be['label']}", **st)
axs[0].set_ylabel("Temperature [K]")
axs[0].set_title("Temperature evolution")
axs[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 4, fontsize = 10)

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
axs[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 4, fontsize=10)

# 3) Stiffness per backend
for be in BACKENDS:
    st = be["style"]
    axs[2].semilogy(time, store[be["key"]]["kappa"],
                    label=f"{be['label']}", **st)
axs[2].set_ylabel(r"Stiffness, $\tau_{max}/\tau_{min}$")
axs[2].set_xlabel("Time [s]")
axs[2].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 2)

# 4) Stiffness ratios vs baseline
#--new code Non-normality
baseline = store[BASELINE_KEY]["kappa"]
for be in BACKENDS:
    st = be["style"]
    #ratio = store[be["key"]]["kappa"] / baseline
    ratio = store[be["key"]]["departure"]
    #axs[3].plot(time, ratio, label=f"{be['label']} / {BASELINE_KEY}", **st)
    axs[3].plot(time, ratio, label=f"{be['label']}", **st)
axs[3].set_ylabel("Nonnormality, AA^T-A^tA")
axs[3].set_xlabel("Time [s]")
axs[3].legend()
plt.show()

