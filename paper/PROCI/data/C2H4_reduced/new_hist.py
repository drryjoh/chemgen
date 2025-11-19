#!python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- Load ---
time = np.load("data/time_s.npy") * 1e6
formulations       = ["cgc", "cgt"]
formulations_label = ["Conservative", "Temperature"]

# --- Style ---
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.labelsize": 16,
    "font.size": 16,
    "legend.fontsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "axes.titlesize": 20,
})

# Okabe–Ito palette
okabe_ito = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#E69F00", "#000000"]
colors = okabe_ito[:2]
markers = ["s", "d"]

# --- Times of interest ---
times_I_want = [0.25, .5, 1.5, 2.5, 4, 8]  # [µs]
t_idx_list = [int(np.argmin(np.abs(time - t))) for t in times_I_want]

# Determine number of modes from one sample file
sample_idx = t_idx_list[0]
modes = len(np.load(f"data/eigs/alphas_{formulations[0]}_{sample_idx}.npy"))
x = np.arange(modes)

# --- Figure with three panels ---
fig, axes = plt.subplots(3, 2, figsize=(10, 10), sharex=True)
axes = axes.flatten()

#handles = [Patch(color=c, label=l) for c, l in zip(colors, formulations_label)]
for ax, t_req, idx_t in zip(axes, times_I_want, t_idx_list):
    for f_idx, (lab, col, mk) in enumerate(zip(formulations_label, colors, markers)):
        y = np.load(f"data/eigs/alphas_{formulations[f_idx]}_{idx_t}.npy").astype(float)
        xvals = np.load(f"data/eigs/svs_{formulations[f_idx]}_{idx_t}.npy").astype(float)
        y_plot = np.where(y > 0.0, y, np.nan)
        ax.loglog(xvals, y_plot, linestyle="none", marker=mk, markersize=10.0, mew=2.0,
                  mec=col, mfc='none', label=lab, alpha=0.9)

    ax.text(0.98, 0.02, f"{t_req:g} µs",
        transform=ax.transAxes,
        ha="right", va="bottom", fontsize=20)
    ax.grid(True, which="both", ls=":", lw=0.6, alpha=0.6)

# Legend at the top center of the first subplot
axes[0].legend(frameon=True, loc="upper left",ncol=1)
axes[4].set_xlabel("$\sigma_i$")
axes[5].set_xlabel("$\sigma_i$")
axes[0].set_ylabel("$\\alpha_i$")
axes[2].set_ylabel("$\\alpha_i$")
axes[4].set_ylabel("$\\alpha_i$")

fig.tight_layout()
plt.savefig("svd_alpha_loglog_points.png", dpi=300)
plt.show()
