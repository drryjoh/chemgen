#!python3
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Load data ----------------
# time_s.npy: time in seconds; converted to microseconds for plotting
time =  np.load("data/time_s.npy") * 1e6            # [µs]
T_ct =  np.load("data/T_fine.npy")                  # [K]
time_ct = np.load("data/time_fine.npy") * 1e6       # [µs]

formulations       = ["cgc","cgt"]
formulations_label = ["conservative", "temperature",]
markers            = ["s","+"]

T_preds = [np.load(f"data/{f}_T_pred.npy") for f in formulations]
kappas  = [np.load(f"data/{f}_kappa.npy")   for f in formulations]

# ---------------- Style (sizes only; keep your existing font family) ----------------
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


# Okabe–Ito colorblind-friendly palette
okabe_ito = [
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#E69F00", "#000000"
]
colors = okabe_ito[:2]  # one color per formulation, reused across T and κ

# ---------------- Figure ----------------
fig, ax1 = plt.subplots(figsize=(6.0, 4.0))  # single-column friendly



temp_handles = []
for lab, m, col, Tpred in zip(formulations_label, markers, colors, T_preds):
    h = ax1.plot(
        time, Tpred,
        linestyle="none", marker=m, ms=5.5, mew=1.1,
        markerfacecolor="none", markeredgecolor=col,
        label=f"Temperature, {lab}"
    )[0]
    temp_handles.append(h)
# Temperature (left y)
temp_ref, = ax1.plot(
    time_ct, T_ct, "-", lw=2.2, color="#000000",
    label="Temperature (cantera)"
)

ax1.set_xlabel(r"time [$\mu$s]")
ax1.set_ylabel("Temperature [K]")
ax1.set_xlim(0, 15)
ax1.grid(True, which="major", alpha=0.25)

# Stiffness (right y, log scale) via twinx
ax2 = ax1.twinx()
ax2.set_yscale("log")

stiff_handles = []
for lab, col, kap in zip(formulations_label, colors, kappas):
    h = ax2.plot(
        time, kap,
        linestyle="--", lw=2.0, color=col,
        label=f"Stiffness, {lab}"
    )[0]
    stiff_handles.append(h)

ax2.set_ylabel(r"Stiffness Ratio $\left(\frac{\lambda_{max}}{\lambda_{min}}\right)_{temp}/\left(\frac{\lambda_{max}}{\lambda_{min}}\right)_{cons}$")

# Keep axes neutral (no custom axis colors)
ax1.tick_params(axis='both', which='both', direction='out')
ax2.tick_params(axis='y', which='both', direction='out')

# Legend: two columns, top-right
handles = [temp_ref] + temp_handles + stiff_handles
labels  = [h.get_label() for h in handles]
ax1.legend(handles, labels, ncol=1, loc="right", frameon=True)

fig.tight_layout()

# Save
plt.savefig("temperature_stiffness_combined.png", dpi=300, bbox_inches="tight")
plt.savefig("temperature_stiffness_combined.pdf", bbox_inches="tight")
plt.show()
