#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

try:
    plt.style.use('seaborn-colorblind')
except OSError:
    plt.style.use('seaborn-v0_8-colorblind')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
mechs = ["OConnaire", "burke", "gri30", "sandiego", "FFCM2_model"]
names = ["Ó Connaire", "Burke", "GRI v3.0", "UCSD", "FFCM2"]
mean_loc = -2.0  # for degree 7

deg = 7
fig, ax = plt.subplots(figsize=(8, 6), sharey=True)
means = []
stds = []

for i, mech in enumerate(mechs):
    path = f"{deg}/l2_{mech}.csv"
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    errors = data[:, 1]
    log_errors = np.log10(errors)

    sns.histplot(
        log_errors,
        ax=ax,
        bins=300,
        stat='density',
        color=colors[i],
        alpha=0.5,
        edgecolor='black',
        label=names[i],
        kde=False  # disable internal kde
        )

# Add KDE manually
    sns.kdeplot(
        log_errors,
        ax=ax,
        color='black',
        linewidth=0.75,
        label=None
        )

    mean_val = np.mean(log_errors)
    std_val = np.std(log_errors)
    means.append(mean_val)
    stds.append(std_val)
    ax.axvline(mean_val, color=colors[i], linestyle='--', linewidth=1)

#ax.set_ylim([0, 350])
txtloc = np.linspace(0.57, 0.9, len(means))
for r, mean_val in enumerate(means):
    ax.text(mean_loc, ax.get_ylim()[1]*(txtloc[r]+0.05), 
            f"{names[r]}: $\mu = {10**mean_val:3.2e}$\n$\sigma = {10**(mean_val + stds[r]):3.2e}$",
            color="k", rotation=0, verticalalignment='top', fontsize=10)

xticks = ax.get_xticks()
integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]
ax.set_xticks(integer_ticks)
ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in integer_ticks])
#ax.set_title(f"$n_p={deg}$", fontsize=14)
ax.set_xlim([-9, 1])
ax.set_ylabel("Probability Density",fontsize=12)
ax.set_xlabel(r"$\epsilon$", fontsize=10)

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles[:5], labels[:5], loc='upper center', ncol=5, fontsize=10, title="Mechanism")

plt.tight_layout(rect=[0, 0.05, 1, 0.90])
plt.savefig("hist.png", dpi=300)
