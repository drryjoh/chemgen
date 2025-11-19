#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm

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

# Global range
vals = [np.load(f"phi_{phi}_p{p}_stiffness.npy") for phi in phis for p in pressures_atm]
all_vals = np.concatenate([v.ravel() for v in vals])

vmin = 1.0              # explicit floor
vmax = 3*10**5            # keep your upper bound

# Grids
nT, nTau = vals[0].shape
T_grid  = np.linspace(1800, 2800, nT)[::-1]
tau_grid = np.linspace(0, 1, nTau)

# Colormap tuned for high-end contrast; mark <1 with 'under' color
cmap = mpl.cm.get_cmap("inferno").copy()
cmap.set_under("#3B4CC0")   # darker, clearer visual block
cmap.set_bad("#F5F5F5")     # NaNs, if any

norm = LogNorm(vmin=vmin, vmax=vmax)

fig, axes = plt.subplots(len(phis), len(pressures_atm), sharex=True, sharey=True)

for i, phi in enumerate(phis):
    for j, p_atm in enumerate(pressures_atm):
        y = np.load(f"phi_{phi}_p{p_atm}_stiffness.npy")  # do NOT clamp; let <1 show as 'under'
        print(np.shape(y))
        im = axes[i, j].pcolormesh(
            tau_grid, T_grid, y,
            shading='nearest', cmap=cmap, norm=norm
        )
        axes[i, j].set_title(rf"$\phi={phi}$,  $p_o={p_atm}$ atm")
        axes[i,j].set_yticks([1800,2800])
        if j == 0:
            axes[i, j].set_ylabel("$T_{o}$ [K]")
        if i == len(phis) - 1:
            axes[i, j].set_xlabel(r"$\tau$")

# Horizontal colorbar on top with log ticks; show under-range triangle
cax = fig.add_axes([0.14, 0.82, 0.72, 0.03])
cbar = fig.colorbar(im, cax=cax, orientation='horizontal', extend='min')
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')
cbar.set_label(r"Stiffness Ratio $\frac{\left(|\Re(\lambda)|_{\max} / |\Re(\lambda)|_{\min}\right)_{\text{temp}}}{\left(|\Re(\lambda)|_{\max} / |\Re(\lambda)|_{\min}\right)_{\text{cons}}}$", fontsize=16, labelpad=12)

# Log tick locator/formatter for clarity; add a few subticks
cbar.locator = mpl.ticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0))
cbar.formatter = mpl.ticker.LogFormatterMathtext(base=10)
cbar.update_ticks()

# Spacing
fig.subplots_adjust(hspace=0.35, wspace=0.25, top=0.75)
plt.savefig("stiffness_parameter_space.png",dpi = 300)
plt.show()
