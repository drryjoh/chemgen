#!python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.patches import Patch

# --- Load ---
time   = np.load("data/time_s.npy") * 1e6
T_ct   = np.load("data/T_fine.npy")
time_ct= np.load("data/time_fine.npy") * 1e6

formulations       = ["cgc","cgt"]
formulations_label = ["conservative", "temperature"]
T_preds = [np.abs(np.load(f"data/{f}_shur_alignment.npy")) for f in formulations]

# --- Style (sizes only) ---
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

# --- 3D grouped bars ---
times_I_want = [0.0, 2.5, 5.0]  # [µs]
t_idx_list   = [int(np.argmin(np.abs(time - t))) for t in times_I_want]

modes = T_preds[0].shape[1]  # e.g., 30
x_centers = np.arange(modes)  # mode indices
y_rows = np.arange(len(times_I_want))  # 0..T-1 for bar placement along y

# bar sizes along x and y (thickness on the floor)
dx_form = 0.38     # x-width of each bar
dy_time = 0.6      # y-depth common to all bars at a given time row

fig = plt.figure(figsize=(7.5, 5.0))
ax = fig.add_subplot(111, projection='3d')

for f_idx, (lab, col, Tpred) in enumerate(zip(formulations_label, colors, T_preds)):
    # offset each formulation left/right within each mode
    x_pos = x_centers + (f_idx - 0.5) * dx_form

    for r, idx_t in enumerate(t_idx_list):
        z0  = np.zeros_like(x_centers, dtype=float)
        dz  = np.load(f"data/eigs/alphas_{formulations[f_idx]}_{idx_t}.npy")

        # Broadcast into vectors for bar3d
        xs = x_pos.astype(float)
        ys = np.full_like(x_centers, r, dtype=float)

        ax.bar3d(xs, ys, z0, dx_form, dy_time, dz, color=col, shade=True, alpha=1.0, linewidth=0)

# Axes and ticks
ax.set_xlabel("Mode index")
ax.set_ylabel("Time [µs]")
ax.set_zlabel("Contribution")

ax.set_xticks(np.arange(0, modes, max(1, modes // 10)))
ax.set_yticks(y_rows)
ax.set_yticklabels([f"{t:g}" for t in times_I_want])

# Optional: tighter z limit if contributions are normalized
# ax.set_zlim(0, 1.0)

# Legend
handles = [Patch(color=c, label=l) for c, l in zip(colors, formulations_label)]
ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.05, 1.0))

plt.savefig("shur_3d_bars.png", dpi=300)
plt.show()
