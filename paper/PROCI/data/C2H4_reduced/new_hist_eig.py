#!python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- Load ---
time = np.load("data/time_s.npy") * 1e6
formulations       = ["cgc", "cgt"]
formulations_label = ["conservative", "temperature"]

# --- Style ---
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.labelsize": 12,
    "font.size": 11,
    "legend.fontsize": 10,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.titlesize": 13,
})

# Okabe–Ito palette
okabe_ito = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#E69F00", "#000000"]
colors = okabe_ito[:2]
markers = ["s", "o"]

# --- Times of interest ---
times_I_want = [0.0, 2.5, 5.0]  # [µs]
t_idx_list = [int(np.argmin(np.abs(time - t))) for t in times_I_want]

# Determine number of modes from one sample file
sample_idx = t_idx_list[0]
modes = len(np.load(f"data/eigs/alphas_{formulations[0]}_{sample_idx}.npy"))
x = np.arange(modes)

# --- Figure with three panels ---
fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), sharey=True)

handles = [Patch(color=c, label=l) for c, l in zip(colors, formulations_label)]
for ax, t_req, idx_t in zip(axes, times_I_want, t_idx_list):
    for f_idx, (lab, col, mk) in enumerate(zip(formulations_label, colors, markers)):
        y = np.load(f"data/eigs/alignment_{formulations[f_idx]}_{idx_t}.npy").astype(float)
        xvals = np.real(np.load(f"data/eigs/eigs_{formulations[f_idx]}_{idx_t}.npy"))
        y_plot = np.where(y > 0.0, y, np.nan)
        ax.loglog(xvals, y_plot, linestyle="none", marker=mk, markersize=4.0,
                  color=col, label=lab, alpha=0.9)

    ax.set_title(f"{t_req:g} µs")
    ax.set_xlabel("Singular value")
    ax.grid(True, which="both", ls=":", lw=0.6, alpha=0.6)

axes[0].set_ylabel("Contribution")

# Legend at the top center of the first subplot
axes[0].legend(handles=handles, frameon=True, loc="upper left",ncol=1)
axes[1].legend(handles=handles, frameon=True, loc="upper left",ncol=1)
axes[2].legend(handles=handles, frameon=True, loc="upper left",ncol=1)

fig.tight_layout()
plt.savefig("svd_alpha_loglog_points.png", dpi=300)
plt.show()
