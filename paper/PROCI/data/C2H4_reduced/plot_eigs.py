#!python3
import numpy as np
import matplotlib.pyplot as plt
# ---------------- Style (sizes only; keep your existing font family) ----------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "axes.labelsize": 14,
    "font.size": 20,
    "legend.fontsize": 20,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.titlesize": 20,
})
# --- Load both datasets ---
n_steps = 200
dt = 2e-7
time = np.linspace(0.0, n_steps * dt, n_steps)
time_us = time * 1e6

times_to_plot  = [1.03, 2.62, 8]

times_I_want = []
for tplt in times_to_plot:
    times_I_want.append(np.argmin(np.abs(tplt-time_us)))

# Colormap by index
# Markers for the two formulations
markers = {"cgc": "s", "cgt": "d"}   # square, triangle
fig, axes = plt.subplots(3,1,figsize=(10,6),sharex=True, sharey=True)
for i, k in enumerate(times_I_want):
    data_cgc = np.load(f"data/eigs/eigs_cgc_{k}.npy")
    data_cgt = np.load(f"data/eigs/eigs_cgt_{k}.npy")
    lams_cgc = data_cgc[np.abs(np.real(data_cgc))>1]
    lams_cgt = data_cgt[np.abs(np.real(data_cgt))>1]
    # --- Temperature formulation (cgc) ---
    axes[i].scatter(
        1/np.abs(np.real(lams_cgc)), np.abs(np.imag(lams_cgc)),
        s=70,
        facecolors="none",
        edgecolors="r",
        linewidths=2,
        marker=markers["cgc"],
        label="conservative",
    )

    axes[i].scatter(
        1/np.abs(np.real(lams_cgt)), np.abs(np.imag(lams_cgt)),
        s=70,
        facecolors="none",
        edgecolors="g",
        linewidths=2,
        marker=markers["cgt"],
        label="temperature",
    )

    axes[i].text(
        0.98, 0.98,
        rf"$t = {times_to_plot[i]:.2f}\,\mu s$",
        transform=axes[i].transAxes,
        ha="right", va="top",
        fontsize=10,
    )

    # Symmetric log scaling
for ax in axes:
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1e-2)

    ax.set_xlabel(r"$1/|\Re(\lambda)|$")
    ax.set_ylabel(r"$\Im(\lambda)$")
    ax.grid(True, alpha=0.3)

# Single legend: remove duplicates
handles, labels = ax.get_legend_handles_labels()
unique = dict(zip(labels, handles))


fig.tight_layout()
plt.show()
