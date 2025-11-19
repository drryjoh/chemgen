#!python3
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Load data ----------------
time = np.load("data/time_s.npy") * 1e6            # [µs]
lin_iters  = [np.load(f"data/{f}_lin_iter.npy")  for f in ["cgc","cgt"]]
newt_iters = [np.load(f"data/{f}_newt_iter.npy") for f in ["cgc","cgt"]]

labels = ["conservative", "temperature"]

# ---------------- Style ----------------
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

# Okabe–Ito (first two colors)
colors = ["#0072B2", "#D55E00"]  # blue, vermillion

# ---------------- Figure ----------------
fig, ax1 = plt.subplots(figsize=(6.0, 4.0))

# Newton iterations (solid lines)
handles1 = []
for lab, col, y in zip(labels, colors, newt_iters):
    h, = ax1.plot(time, y, linestyle="-", linewidth=2.0, color=col,
                  label=f"Newton iterations, {lab}")
    handles1.append(h)

ax1.set_xlabel(r"time [$\mu$s]")
ax1.set_ylabel("Total Newton iterations per time step")
ax1.set_xlim(0, 15)
ax1.set_ylim(0, max(max(v) for v in newt_iters)*1.1)
ax1.grid(True, which="major", alpha=0.25)

# Linear iterations / stiffness (dashed lines), same colors
ax2 = ax1.twinx()
handles2 = []
for lab, col, y in zip(labels, colors, lin_iters):
    h, = ax2.plot(time, y, linestyle="--", linewidth=2.0, color=col,
                  label=f"Linear iterations, {lab}")
    handles2.append(h)
idx_time =  np.argmin(np.abs(time-15))
print(f"Total number fo linear iterations for {labels[0]}: {np.sum(lin_iters[0])}")
print(f"Total number fo linear iterations for {labels[1]}: {np.sum(lin_iters[1])}")
ax2.set_ylabel("Total linear iterations per time step")
ax2.set_ylim(0, max(max(v) for v in lin_iters)*1.1)

# Legend (combine)
handles = handles1 + handles2
labels_leg = [h.get_label() for h in handles]
ax1.legend(handles, labels_leg, ncol=1, loc="upper right", frameon=True)

fig.tight_layout()
plt.savefig("iterations_combined.png", dpi=300, bbox_inches="tight")
plt.savefig("iterations_combined.pdf", bbox_inches="tight")
plt.show()
