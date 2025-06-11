#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Plot configuration
colors = ["purple", "blue", "orange", "red", "green"]
mechs = ["OConnaire", "burke", "gri30", "sandiego", "FFCM2_model"]
names = ["Ó Connaire", "Burke", "GRI v3.0", "UCSD", "FFCM2"]
degrees = [4, 7, 15]  # thermodynamic polynomial degrees
mean_loc = [-8.5,-2,-2]
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

for j, deg in enumerate(degrees):
    ax = axes[j]
    means = []
    stds = []
    for i, mech in enumerate(mechs):
        # Construct file path from directory
        path = f"{deg}/l2_{mech}.csv"
        data = np.loadtxt(path, delimiter=",", skiprows=1)
        errors = data[:, 1]
        log_errors = np.log10(errors)
        
        # Plot histogram
        sns.histplot(
            log_errors,
            ax=ax,
            bins=300,
            color=colors[i],
            alpha=0.5,
            edgecolor='black',
            label=names[i],
            kde=True
        )

        # Plot mean line
        mean_val = np.mean(log_errors)
        std_val = np.std(log_errors)
        means.append(mean_val)
        stds.append(std_val)
        ax.axvline(mean_val, color=colors[i], linestyle='--', linewidth=1)
    ax.set_ylim([0,350])
    txtloc = np.linspace(0.6,0.9,len(means))
    for r, mean_val in enumerate(means):
        ax.text(mean_loc[j], ax.get_ylim()[1]*(txtloc[r]+0.05), f"{names[r]}: $\mu = $ {(10**mean_val):3.2e} \n $\sigma =$ {(10**(mean_val+stds[r])):3.2e}",
        color=colors[r], rotation=0, verticalalignment='top', fontsize=8)

    # X-axis ticks and labels
    xticks = ax.get_xticks()
    integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]
    ax.set_xticks(integer_ticks)
    ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in integer_ticks])
    ax.set_title(f"$n_p={deg}$", fontsize=14)
    ax.set_xlim([-9,1])

    if j == 0:
        ax.set_ylabel("Distribution", fontsize=12)
    ax.set_xlabel(r"$\epsilon$", fontsize=10)

# Extract handles/labels from first subplot
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles[:5], labels[:5], loc='upper center', ncol=5, fontsize=10, title="Mechanism")

#fig.suptitle("Log-Scaled Error Distributions by Mechanism and Polynomial Degree", fontsize=16)
plt.tight_layout(rect=[0, 0.05, 1, 0.90])
plt.savefig("combined_hist.png", dpi=300)
plt.show()
