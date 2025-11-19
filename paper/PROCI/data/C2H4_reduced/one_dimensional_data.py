#!python3
#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.colors as mcolors  # kept only if you later reuse styles
import cantera as ct

# ChemGen backends
import chemgen_conservative as cgc
import chemgen_temperature as cgt

from numerics_refactored import backwards_euler, sdirk2, full_jacobian
from metrics_refactored import (
    inv_largest_tau, inv_smallest_tau, eigs, svs, alphas, stiffness,
    departure, schur_alignment, liao_metric, contribution_to_b, modal_contribution, svd_alignment, PALETTE
)

# ------------------------------- helpers --------------------------------
DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)

def save_npy(name: str, arr):
    """Save to data/<name>.npy."""
    np.save(os.path.join(DATA_DIR, f"{name}.npy"), arr)

def save_npz(name: str, **arrays):
    """Save to data/<name>.npz (compressed)."""
    np.savez_compressed(os.path.join(DATA_DIR, f"{name}.npz"), **arrays)

def stack_or_object(seq_of_arrays):
    """
    Stack per-step arrays to 2D if shapes are consistent.
    Otherwise return object array for safe serialization.
    """
    try:
        return np.stack(seq_of_arrays, axis=0)
    except Exception:
        return np.array(seq_of_arrays, dtype=object)

def trim_outliers_log(vals, n_iter=12):
    """Iteratively remove the element farthest from mean(log(vals))."""
    v = np.array(vals, dtype=float)
    if v.size == 0:
        return v
    w = v.copy()
    for _ in range(min(n_iter, max(0, w.size - 1))):
        m = np.mean(np.log(w))
        idx = np.argmax(np.abs(np.log(w) - m))
        w = np.delete(w, idx)
    return w

def clustering_scores_log(x, alpha=0.3):
    x = np.asarray(x, dtype=float)
    eps = 1e-30
    x_log = np.log10(np.abs(x) + eps)
    x_log = x_log[np.isfinite(x_log)]
    if x_log.size < 2:
        return np.nan, np.nan, np.nan
    x_log.sort()
    d = np.empty_like(x_log)
    d[0]  = x_log[1] - x_log[0]
    d[-1] = x_log[-1] - x_log[-2]
    d[1:-1] = np.minimum(x_log[1:-1] - x_log[:-2], x_log[2:] - x_log[1:-1])
    mean_d = np.mean(d)
    cv_nn  = np.std(d) / mean_d
    norm   = (x_log[-1] - x_log[0]) / (len(x_log) - 1)
    R_nn   = mean_d / norm
    F_eps  = np.mean(d < alpha * norm)
    return R_nn, cv_nn, F_eps

# ---------------------------- backend registry ---------------------------
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
for i, be in enumerate(BACKENDS):
    be["style"] = dict(color=PALETTE[i % len(PALETTE)],
                       ls="--" if i % 2 else "-.",
                       lw=2)
BASELINE_KEY = "cgc" if any(b["key"] == "cgc" for b in BACKENDS) else BACKENDS[0]["key"]

# ------------------------------ Cantera set-up ---------------------------
gas = ct.Solution("mechanism.yaml")
gas.set_equivalence_ratio(phi=1.0, fuel="C2H4", oxidizer="O2:1.0, N2:3.76")
gas.TP = 2500.0, 101325.0
reactor = ct.IdealGasReactor(gas)
net = ct.ReactorNet([reactor])

n_steps = 400
dt = 1e-7
time = np.linspace(0.0, n_steps * dt, n_steps)
time_us = time * 1e6

#fine data
# Build a fresh gas with the SAME mechanism and initial state
gas_fine = ct.Solution("mechanism.yaml")
gas_fine.set_equivalence_ratio(phi=1.0, fuel="C2H4", oxidizer="O2:1.0, N2:3.76")
gas_fine.TP = 2500.0, 101325.0   # match your coarse run's initial state (TPX if needed)

# Independent reactor/network
reactor_fine = ct.IdealGasReactor(gas_fine)
net_fine = ct.ReactorNet([reactor_fine])

# Time grid (include t=0) and storage
time_fine = np.linspace(0.0, n_steps * dt, n_steps * 10 + 1)
T_fine = [reactor_fine.T]  # initial value at t=0

# Advance with ABSOLUTE times using the fine network
for t in time_fine[1:]:
    net_fine.advance(t)
    T_fine.append(reactor_fine.T)

# Save
save_npy("time_fine", time_fine)
save_npy("T_fine", np.array(T_fine))

# Reference trajectory
T_phys = []

# Per-backend storage
store = {
    be["key"]: dict(
        T_pred=[], inv_tau_min=[], inv_tau_max=[], kappa=[],
        lin_iter=[], newt_iter=[], eigs=[], svs=[], departure=[], alphas=[], shur_alignment=[]
    ) for be in BACKENDS
}

# --------------------------- time integration loop -----------------------
for time_indx, t in enumerate(time):
    T = reactor.T
    C = reactor.thermo.concentrations

    # Advance all backends explicitly on same (C,T)
    for be in BACKENDS:
        # y, n_lin, n_newt = sdirk2(C, T, dt, be, gmres_method="numpy")
        y, n_lin, n_newt = backwards_euler(C, T, dt, be, gmres_tolerance=1e-8, gmres_method="numpy", linear_solver_verbose = True, step = time_indx, formulation = be["key"])
        T_hat = be["temp_from_state"](y)
        s = store[be["key"]]
        s["T_pred"].append(T_hat)
        s["lin_iter"].append(n_lin)
        s["newt_iter"].append(n_newt)

    # Advance Cantera to physical time t
    net.advance(t)
    T_phys.append(reactor.T)

    # Jacobians at updated physical state
    T_now = reactor.T
    C_now = reactor.thermo.concentrations
    for be in BACKENDS:
        mod = be["mod"]
        J = np.array(mod.source_jacobian(C_now, T_now))
        res = mod.source(C_now, T_now)
        s = store[be["key"]]
        s["inv_tau_min"].append(inv_largest_tau(J))
        s["inv_tau_max"].append(inv_smallest_tau(J))
        eigs_t = eigs(J)
        s["eigs"].append(eigs_t)
        svs_a = svs(full_jacobian(J, dt))
        s["svs"].append(svs_a)
        alphas_t = svd_alignment(full_jacobian(J, dt/100), np.array(res))
        s["shur_alignment"].append(alphas_t)#modal_contribution(full_jacobian(J, dt/100), -np.array(res), 1, 1)
        eig_align = contribution_to_b(full_jacobian(J, dt), np.array(res))
        np.save(f"data/eigs/eigs_{be['key']}_{time_indx}.npy",eigs_t)
        np.save(f"data/eigs/svs_{be['key']}_{time_indx}.npy",svs_a)
        np.save(f"data/eigs/alphas_{be['key']}_{time_indx}.npy",alphas_t)
        np.save(f"data/eigs/alignment_{be['key']}_{time_indx}.npy",eig_align)
        s["alphas"].append(alphas_t)

        np.save(f"data/eigs/eigs_{be['key']}_{time_indx}.npy",eigs_t)
        s["kappa"].append(stiffness(J))
        s["departure"].append(departure(J))

# ------------------------------- to arrays -------------------------------
T_phys = np.asarray(T_phys)
save_npy("time_s", time)
save_npy("time_us", time_us)
save_npy("T_phys", T_phys)

for be in BACKENDS:
    key = be["key"]
    s = store[key]

    # 1D series
    T_pred      = np.asarray(s["T_pred"])
    inv_tau_min = np.asarray(s["inv_tau_min"])
    inv_tau_max = np.asarray(s["inv_tau_max"])
    kappa       = np.asarray(s["kappa"])
    lin_iter    = np.asarray(s["lin_iter"])
    newt_iter   = np.asarray(s["newt_iter"])
    departure_  = np.asarray(s["departure"])
    shur_alignment  = np.asarray(s["shur_alignment"])

    save_npy(f"{key}_T_pred", T_pred)
    save_npy(f"{key}_inv_tau_min", inv_tau_min)
    save_npy(f"{key}_inv_tau_max", inv_tau_max)
    save_npy(f"{key}_kappa", kappa)
    save_npy(f"{key}_lin_iter", lin_iter)
    save_npy(f"{key}_newt_iter", newt_iter)
    save_npy(f"{key}_departure", departure_)
    save_npy(f"{key}_shur_alignment", shur_alignment)

    # Per-step spectra (may be 1) constant-length -> stackable, or 2) variable-length)
    eigs_all  = stack_or_object(s["eigs"])
    svs_all   = stack_or_object(s["svs"])
    alphas_all= stack_or_object(s["alphas"])

    # Save as compressed .npz to keep keys and support object arrays if needed
    save_npz(f"{key}_eigs",  eigs=eigs_all)
    save_npz(f"{key}_svs",   svs=svs_all)
    save_npz(f"{key}_alphas",alphas=alphas_all)

    # Optional: clustering metrics from trimmed singular values
    R_list, CV_list, F_list = [], [], []
    for sv_vec in s["svs"]:
        sv_trim = trim_outliers_log(sv_vec, n_iter=12)
        R, CV, F = clustering_scores_log(sv_trim, alpha=0.3)
        R_list.append(R); CV_list.append(CV); F_list.append(F)
    save_npy(f"{key}_cluster_RNN", np.asarray(R_list))
    save_npy(f"{key}_cluster_CVNN", np.asarray(CV_list))
    save_npy(f"{key}_cluster_Feps", np.asarray(F_list))

# Baseline reference if needed downstream
save_npy("baseline_key", np.array(Baseline := BASELINE_KEY))

print("Done. Arrays saved under ./data/")
