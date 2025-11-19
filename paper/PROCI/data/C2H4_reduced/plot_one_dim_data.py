#!python3
#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ---- config ----
DATA_DIR = "data"
FIG_DIR = "figs"
os.makedirs(FIG_DIR, exist_ok=True)

# Backends and styles (keys must match what you saved)
BACKENDS = [
    {"key": "cgc", "label": "Conservative"},
    {"key": "cgt", "label": "Temperature"},
]
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

for i, be in enumerate(BACKENDS):
    be["style"] = dict(color=PALETTE[i % len(PALETTE)],
                       ls="--" if i % 2 else "-.",
                       lw=2)

# ---- helpers ----
def load(path):
    return np.load(path) if os.path.exists(path) else None

def draw_density_ellipse_curve(ax, xdata, ydata, color="red", lw=1.5, scale=2.0, npts=256, **plot_kwargs):
    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)
    mask = (xdata > 0) & (ydata > 0) & np.isfinite(xdata) & np.isfinite(ydata)
    if mask.sum() < 3:
        return None
    X = np.log10(xdata[mask])
    Y = np.log10(ydata[mask])
    mean = np.array([X.mean(), Y.mean()])
    cov = np.cov(X, Y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    rx, ry = scale * np.sqrt(np.maximum(vals, 0.0))
    t = np.linspace(0, 2*np.pi, npts)
    R = vecs @ np.diag([rx, ry])
    pts_log = mean[:, None] + R @ np.vstack((np.cos(t), np.sin(t)))
    xl = 10**pts_log[0, :]
    yl = 10**pts_log[1, :]
    return ax.plot(xl, yl, color=color, lw=lw, **plot_kwargs)

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

# ---- load shared series ----
# Prefer a shared time_us.npy if present; otherwise use one backend’s time
time_us = load(os.path.join(DATA_DIR, "time_us.npy"))
if time_us is None:
    for be in BACKENDS:
        t_try = load(os.path.join(DATA_DIR, f"time_{be['key']}.npy"))
        if t_try is not None:
            time_us = t_try
            break
time = load(os.path.join(DATA_DIR, "time.npy"))
if time is None and time_us is not None:
    time = time_us / 1e6

T_phys = load(os.path.join(DATA_DIR, "T_phys.npy"))

# ---- Figure 1: singular values vs time + clustering metrics ----
# expects per-backend svs as a list/array of arrays: svs_{key}.npy
fig, ax = plt.subplots(2, 2, figsize=(15, 10), sharey=False)
axs = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]

if time_us is not None:
    for be in BACKENDS:
        svs_all = load(os.path.join(DATA_DIR, f"svs_{be['key']}.npy"))
        if svs_all is None:
            continue
        st = be["style"]
        marker = "o"
        # svs_all is expected shape [nt] list-like; handle both list-of-arrays or object array
        for i, t_us in enumerate(time_us):
            if i >= len(svs_all):
                break
            svs_t = np.array(svs_all[i]).astype(float).ravel()
            svs_t = svs_t[np.isfinite(svs_t) & (svs_t > 0)]
            if svs_t.size == 0:
                continue
            # optional windowing as in original
            for _ in range(min(12, svs_t.size-1)):
                mean_log = np.mean(np.log(svs_t))
                idx_rem = np.argmax(np.abs(np.log(svs_t) - mean_log))
                svs_t = np.delete(svs_t, idx_rem)
                if svs_t.size < 2:
                    break
            if svs_t.size == 0:
                continue
            axs[0].semilogy(svs_t*0 + t_us, svs_t, marker,
                            mfc=mcolors.to_rgba(st["color"], alpha=0.7),
                            mec='black', markersize=4,
                            label=be["label"] if i == 0 else None)
            cluster = clustering_scores_log(svs_t)
            axs[2].plot(t_us, cluster["R_NN"], 'o',
                        mfc=mcolors.to_rgba(st["color"], alpha=0.7),
                        mec='black', markersize=4,
                        label=be["label"] if i == 0 else None)
            axs[3].plot(t_us, cluster["CV_NN"], 'o',
                        mfc=mcolors.to_rgba(st["color"], alpha=0.7),
                        mec='black', markersize=4,
                        label=be["label"] if i == 0 else None)

for a in axs:
    a.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
    a.grid(True, which="minor", linestyle="--", linewidth=0.3, alpha=0.5)
axs[2].set_ylabel("R_NN")
axs[3].set_ylabel("CV_NN")
axs[0].set_ylabel("Singular values")
axs[0].set_xlabel(r"Time [$\mu$s]")
axs[2].set_xlabel(r"Time [$\mu$s]")
axs[3].set_xlabel(r"Time [$\mu$s]")
axs[0].legend(ncol=2, frameon=False)
axs[2].legend(ncol=2, frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "svs_clustering.png"), dpi=200, bbox_inches="tight")

# ---- Figure 2: eigs vs alphas with density ellipse ----
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
for be in BACKENDS:
    eigs_all = load(os.path.join(DATA_DIR, f"{be['label']}_svs.npy"))
    alphas_all = load(os.path.join(DATA_DIR, f"{be['label']}_alphas.npy"))
    if eigs_all is None or alphas_all is None:
        continue
    st = be["style"]
    x = np.asarray(eigs_all).ravel()
    y = np.asarray(alphas_all).ravel()
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (np.abs(y) > 0)
    x, y = x[mask], np.abs(y[mask])
    ax.loglog(x, y, 'o',
              mfc=mcolors.to_rgba(st["color"], 0.7), mec='black', markersize=3,
              label=be["label"])
    draw_density_ellipse_curve(ax, x, y, color=st["color"], scale=4.0, lw=2.5)

ax.set_xlabel(r"Singular values, $\sigma_i$")
ax.set_ylabel(r"$b$-alignment, $|\alpha_i|$")
ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
ax.legend(frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "eigs_alignment.png"), dpi=300, bbox_inches="tight")

# ---- Figure 3: temperature, tau, stiffness, nonnormality vs time_us ----
fig, axs = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
if T_phys is not None and time_us is not None:
    axs[0].plot(time_us, T_phys, color="#000000", label="Cantera")
for be in BACKENDS:
    if time_us is None:
        break
    st = be["style"]
    T_pred = load(os.path.join(DATA_DIR, f"T_pred_{be['key']}.npy"))
    if T_pred is not None:
        axs[0].plot(time_us, T_pred, label=be["label"], **st)
axs[0].set_ylabel("Temperature [K]")
axs[0].set_title("Temperature evolution")
axs[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=2, fontsize=10)

for be in BACKENDS:
    if time_us is None:
        continue
    st = be["style"]
    inv_tau_min = load(os.path.join(DATA_DIR, f"inv_tau_min_{be['key']}.npy"))
    inv_tau_max = load(os.path.join(DATA_DIR, f"inv_tau_max_{be['key']}.npy"))
    if inv_tau_min is not None:
        axs[1].semilogy(time_us, inv_tau_min, label=rf"$\tau_{{min}}$ {be['label']}", **st)
    if inv_tau_max is not None:
        axs[1].semilogy(time_us, inv_tau_max, linestyle=":", color=st["color"],
                        label=rf"$\tau_{{max}}$ {be['label']}")
axs[1].set_ylabel(r"$\tau = 1/\lambda$")
axs[1].set_title("Inverse of largest and smallest nonzero eigenvalues")
axs[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=2, fontsize=10)

for be in BACKENDS:
    if time_us is None:
        continue
    st = be["style"]
    kappa = load(os.path.join(DATA_DIR, f"kappa_{be['key']}.npy"))
    if kappa is not None:
        axs[2].semilogy(time_us, kappa, label=be["label"], **st)
axs[2].set_ylabel(r"Stiffness, $\tau_{max}/\tau_{min}$")
axs[2].set_xlabel(r"Time [$\mu$s]")
axs[2].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=2)

for be in BACKENDS:
    if time_us is None:
        continue
    st = be["style"]
    dep = load(os.path.join(DATA_DIR, f"departure_{be['key']}.npy"))
    if dep is not None:
        axs[3].plot(time_us, dep, label=be["label"], **st)
axs[3].axhline(1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.8)
axs[3].set_ylabel("Nonnormality proxy")
axs[3].set_xlabel(r"Time [$\mu$s]")
axs[3].legend(frameon=False)

for ax in axs:
    ax.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
    ax.grid(True, which="minor", linestyle="--", linewidth=0.3, alpha=0.5)

plt.tight_layout()
plt.subplots_adjust(right=0.78)
plt.savefig(os.path.join(FIG_DIR, "stiffness_formulation.png"), dpi=300, bbox_inches="tight")

# ---- Figure 4: iteration metrics ----
fig, axs = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
if time_us is not None:
    for be in BACKENDS:
        st = be["style"]
        lin = load(os.path.join(DATA_DIR, f"linear_iterations_{be['key']}.npy"))
        if lin is not None:
            total = int(np.nansum(lin))
            axs[0].plot(time_us, lin, label=f"{be['label']} total = {total}", **st)
    axs[0].set_ylabel("Linear iters / step")
    axs[0].legend(ncol=2, frameon=False)

    for be in BACKENDS:
        st = be["style"]
        newt = load(os.path.join(DATA_DIR, f"newt_iter_{be['key']}.npy")) \
               or load(os.path.join(DATA_DIR, f"newton_iterations_{be['key']}.npy"))
        if newt is not None:
            axs[1].plot(time_us, newt, label=be["label"], **st)
    axs[1].set_ylabel("Nonlinear iters / step")
    axs[1].legend(ncol=2, frameon=False)

    for be in BACKENDS:
        st = be["style"]
        lin = load(os.path.join(DATA_DIR, f"linear_iterations_{be['key']}.npy"))
        newt = load(os.path.join(DATA_DIR, f"newt_iter_{be['key']}.npy")) \
               or load(os.path.join(DATA_DIR, f"newton_iterations_{be['key']}.npy"))
        if lin is None or newt is None:
            continue
        lin = np.asarray(lin, float)
        newt = np.asarray(newt, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(newt > 0, lin / newt, np.nan)
        axs[2].plot(time_us, ratio, label=be["label"], **st)
    axs[2].set_ylabel("Linear iters / nonlinear step")
    axs[2].set_xlabel(r"Time [$\mu$s]")
    axs[2].legend(ncol=2, frameon=False)

    for ax in axs:
        ax.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
        ax.grid(True, which="minor", linestyle="--", linewidth=0.3, alpha=0.5)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "iteration_metrics.png"), dpi=300, bbox_inches="tight")

print("Saved figures:",
      "figs/svs_clustering.png, figs/eigs_alignment.png, figs/stiffness_formulation.png, figs/iteration_metrics.png")

