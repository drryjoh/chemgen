#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Style
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Garamond", "Palatino Linotype"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "figure.figsize": (7.0, 6.0),
    "axes.titlesize": 16,
    "lines.linewidth": 1.2,
    "axes.grid": False,
})

pressures_atm = [1, 5, 10]
phis = [0.8, 1.0, 1.2]

# Global range (linear)
liner_or_nonlinear = "nonlinear"
vals = [np.load(f"phi_{phi}_p{p}_{liner_or_nonlinear}.npy") for phi in phis for p in pressures_atm]
all_vals = np.concatenate([v.ravel() for v in vals])
print(np.mean(all_vals))
print(np.median(all_vals))
print(np.std(all_vals))

vmin = 1
vmax = 2

# Grids
nT, nTau = vals[0].shape
T_grid  = np.linspace(1800, 2800, nT)[::-1]
tau_grid = np.linspace(0, 1, nTau)

# Colormap: inferno + explicit "over" color for > vmax
cmap = mpl.cm.get_cmap("viridis").copy()
cmap.set_under("red")   # distinct blue for values > vmax
cmap.set_over("black")   # distinct blue for values > vmax
cmap.set_bad("#F5F5F5")

fig, axes = plt.subplots(len(phis), len(pressures_atm), sharex=True, sharey=True)

for i, phi in enumerate(phis):
    for j, p_atm in enumerate(pressures_atm):
        y = np.load(f"phi_{phi}_p{p_atm}_{liner_or_nonlinear}.npy")

        im = axes[i, j].pcolormesh(
            tau_grid, T_grid, y,
            shading='nearest', cmap=cmap,
            vmin=vmin, vmax=vmax
        )

        # Optional: draw contour at the threshold (vmax)
        #axes[i, j].contour(tau_grid, T_grid, y, levels=[vmax],
        #                   colors='k', linewidths=0.8)

        axes[i, j].set_title(rf"$\phi={phi}$, $p_o={p_atm}$ atm")
        axes[i, j].set_yticks([1800, 2800])
        if j == 0:
            axes[i, j].set_ylabel("$T_{o}$ [K]")
        if i == len(phis) - 1:
            axes[i, j].set_xlabel(r"$\tau$")

# Horizontal colorbar on top; use 'extend=max' to show the >vmax triangle
cax = fig.add_axes([0.14, 0.82, 0.72, 0.03])
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')
cbar.set_label(
    r"$n_{\text{nl, temp}}/n_{\text{nl, cons}}$",
    fontsize=16, labelpad=10
)

fig.subplots_adjust(hspace=0.35, wspace=0.25, top=0.75)
plt.savefig("stiffness_parameter_space_nonlinear_threshold.png", dpi=300)
plt.show()
