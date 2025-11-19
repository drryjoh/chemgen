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

T_preds = [np.load(f"data/{f}_shur_alignment.npy") for f in formulations]
for Tpredi in T_preds:
    print(Tpredi)

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
four_colors = okabe_ito[:4]

# ---------------- Figure ----------------
# ---------------- Figure ----------------
fig, ax1 = plt.subplots(figsize=(6.0, 4.0))  # single-column friendly

temp_handles = []
cidx = 0
for lab, m, col, Tpred in zip(formulations_label, markers, colors, T_preds):
    # Lower spectral band V1–V5
    h = ax1.semilogy(
        time, Tpred[:,0],
        linestyle="none", marker=m, ms=5.5, mew=1.1,
        markerfacecolor="none", markeredgecolor=four_colors[cidx],
        label=rf"{lab}: Liao"#sum $\eta_1$ to $\eta_{{5}}$"
    )[0]
    temp_handles.append(h)
    cidx += 1

    ## Upper spectral band \eta_{n-5}–\eta_n
    #h = ax1.plot(
    #    time, Tpred[:,1],
    #    linestyle="none", marker=m, ms=5.5, mew=1.1,
    #    markerfacecolor="none", markeredgecolor=four_colors[cidx],
    #    label=rf"{lab}: $\eta_n$"
    #)[0]
    #temp_handles.append(h)
    #cidx += 1

ax1.set_xlabel(r"time [$\mu$s]")
ax1.set_ylabel("Schur mode alignment")
ax1.set_xlim(0, 15)
ax1.grid(True, which="major", alpha=0.25)

# ---- Legend above plot ----
ax1.legend(
    ncol=2,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.00),
    frameon=True,
)

fig.tight_layout()
plt.savefig("shur.png",dpi=300)
plt.show()
