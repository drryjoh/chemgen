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

from numerics_refactored import backwards_euler, sdirk2, full_jacobian
from metrics_refactored import (inv_largest_tau, inv_smallest_tau, eigs, svs, alphas, stiffness,
                     departure, PALETTE)
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
                         lin_iter=[], newt_iter=[], eigs=[], n_eigs=[], svs=[], departure=[], alphas=[])
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
        n_eigs, a_eigs = eigs(full_jacobian(J, dt*1000))
        s["eigs"].append(a_eigs)
        s["n_eigs"].append(n_eigs)
        s["svs"].append(svs(full_jacobian(J, dt)))
        s["alphas"].append(alphas(J, be["mod"].source(C_now, T_now)))
        s["kappa"].append(stiffness(J))
        s["departure"].append(departure(J))

n_cgt_eigs = []
n_cgt_eigs = []
plt.figure()
for be in BACKENDS:
    s = store[be["key"]]
    a_eigs = s["eigs"]
    print(np.sort(np.abs(np.imag(s["alphas"][10]/(np.real(s["alphas"][10])+1e-16)))))
    n_eigs = s["n_eigs"]
    plt.plot(time, n_eigs)
plt.show()
