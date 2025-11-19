#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def truncate_cmap(cmap, minval=0.0, maxval=1.0, n=256):
    """Return a copy of cmap truncated to [minval, maxval]."""
    new_colors = cmap(np.linspace(minval, maxval, n))
    return mcolors.ListedColormap(new_colors)

# === Load data ===
csvs = np.load("Conservative_svs.npy")
tsvs = np.load("Temperature_svs.npy")
calphas = np.load("Conservative_alphas.npy")
talphas = np.load("Temperature_alphas.npy")

# === Global style ===
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Garamond", "Times New Roman", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11
})

# === Create figure ===
fig, ax = plt.subplots(1, 1, figsize=(20, 8))
axes = [ax]
# use base colormaps
base_cons = plt.cm.Blues
base_temp = plt.cm.YlOrBr

# truncate to match color_idx start at 0.4
cmap_cons = truncate_cmap(base_cons, 0.4, 1.0)
cmap_temp = truncate_cmap(base_temp, 0.4, 1.0)

# now use full [0,1] for indices because the maps are already truncated
color_idx = np.linspace(0.0, 1.0, len(csvs))

# === Plot data ===
for i, (csv, tsv, ca, ta) in enumerate(zip(csvs, tsvs, calphas, talphas)):
    color_c = cmap_cons(color_idx[i])
    color_t = cmap_temp(color_idx[i])
    
    if i==0:
        axes[0].loglog(csv, np.abs(ca), 'o', linestyle='none',
                    markerfacecolor="k", markeredgewidth=0, markersize=3)
        print(np.sort(csv))
        axes[0].loglog(tsv, np.abs(ta), '^', linestyle='none',
                    markerfacecolor="r", markeredgewidth=0, markersize=3)
        print(np.sort(tsv))

    #axes[1].loglog(csv, np.abs(ca), 'o', linestyle='none',
    #               markerfacecolor=color_c, markeredgewidth=0, markersize=3)
    #axes[1].loglog(tsv, np.abs(ta), '^', linestyle='none',
    #               markerfacecolor=color_t, markeredgewidth=0, markersize=3)

# === Zoomed-in limits ===
#xlim_zoom = [2.5e5, 5e9]
#ylim_zoom = [2e-6, 2e3]
#axes[1].set_xlim(xlim_zoom)
#axes[1].set_ylim(ylim_zoom)

## === Rectangle on left plot ===
#rect = Rectangle((xlim_zoom[0], ylim_zoom[0]),
#                 xlim_zoom[1] - xlim_zoom[0],
#                 ylim_zoom[1] - ylim_zoom[0],
#                 linewidth=1.0, edgecolor='black', facecolor='none', linestyle="--")
#axes[0].add_patch(rect)

# === Axis labels and grid ===
for ax in axes:
    ax.set_xlabel(r"Eigen values, $\lambda_i$")
    ax.set_ylabel(r"$\alpha_i$")
    ax.grid(True, which="both", ls=":", lw=0.5)

# === Add two small horizontal colorbars INSIDE axes[0] ===
tau_label = r"$\tau$"
# Match color range to what you're actually plotting (0.4–1.0)
vmin, vmax = 0.4, 1.0
norm = Normalize(vmin=vmin, vmax=vmax)

sm_temp = cm.ScalarMappable(norm=norm, cmap=cmap_temp)
sm_cons = cm.ScalarMappable(norm=norm, cmap=cmap_cons)


# Dimensions in axis fraction coordinates
cb_width = 0.25
cb_height = 0.02
pad = 0.05  # vertical separation

# Conservative (bottom bar)
cax_cons = axes[0].inset_axes([0.07, 0.8, cb_width, cb_height])
cb_cons = plt.colorbar(sm_cons, cax=cax_cons, orientation='horizontal')
cb_cons.set_ticks([0.4, 1])
cb_cons.set_ticklabels(["0", tau_label])
cb_cons.ax.tick_params(labelsize=10)
# Label above bar
cax_cons.text(0.5, 1., "Conservative Formulation", ha='center', va='bottom',
              transform=cax_cons.transAxes, fontsize=11, fontstyle='italic')

# Temperature (top bar)
cax_temp = axes[0].inset_axes([0.07, 0.8 + cb_height + pad, cb_width, cb_height])
cb_temp = plt.colorbar(sm_temp, cax=cax_temp, orientation='horizontal')
cb_temp.set_ticks([0.4, 1])
cb_temp.set_ticklabels(["$0$", tau_label])
cb_temp.ax.tick_params(labelsize=10)
# Label above bar
cax_temp.text(0.5, 1., "Temperature Formulation", ha='center', va='bottom',
              transform=cax_temp.transAxes, fontsize=11, fontstyle='italic')

plt.tight_layout()
plt.savefig("sv_vs_alpha.png",dpi=300)
plt.show()
