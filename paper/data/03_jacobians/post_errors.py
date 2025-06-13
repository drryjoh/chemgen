#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Setup
fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

colors = ["purple", "blue", "orange", "red", "green"]
mechs = ["OConnaire", "burke", "gri30", "sandiego", "FFCM2_model"]
names = ["Ó Connaire", "Burke", "GRI v3.0", "UCSD", "FFCM2"]

L2s = [f"errors/L2_{mech}.npy" for mech in mechs]
L2snei = [f"errors/L2_nei_{mech}.npy" for mech in mechs]

# Left subplot: pointwise error
ax = axes[0]
print("LAPACK Error Summary:")
print(f"{'Mechanism':<12} {'Mean Error':>12} {'Max Error':>12}")
for k, L2 in enumerate(L2s):
    errors = np.load(L2)
    log_errors = np.log10(errors)
    
    mean_val = np.mean(errors)
    max_val = np.max(errors)
    print(f"{names[k]:<12} {mean_val:12.2e} {max_val:12.2e}")

    sns.histplot(
        log_errors,
        ax=ax,
        bins=300,
        color=colors[k],
        alpha=0.5,
        edgecolor='black',
        label=names[k],
        kde=True
    )
    ax.axvline(np.mean(log_errors), color=colors[k], linestyle='--', linewidth=1)

ax.set_xlim([-14,-2])
xticks = ax.get_xticks()
integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]
ax.set_xticks(integer_ticks)
ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in integer_ticks])
ax.set_title("$E_{norm}$ Error")
ax.set_xlabel("$E_{norm}$ Error")
ax.set_ylabel("Frequency")
ax.legend(title="Mechanism")

# Right subplot: neighbor-aware error
ax = axes[1]
print("\nNeimeyer ErrorError Summary:")
print(f"{'Mechanism':<12} {'Mean Error':>12} {'Max Error':>12}")
for k, L2 in enumerate(L2snei):
    errors = np.load(L2)
    log_errors = np.log10(errors)

    mean_val = np.mean(errors)
    max_val = np.min(np.sort(errors)[-5:][::-1])
    print(f"{names[k]:<12} {mean_val:12.2e} {max_val:12.2e}")


    sns.histplot(
        log_errors,
        ax=ax,
        bins=300,
        color=colors[k],
        alpha=0.5,
        edgecolor='black',
        label=names[k],
        kde=True
    )
    ax.axvline(np.mean(log_errors), color=colors[k], linestyle='--', linewidth=1)
ax.set_xlim([-14,-2])
xticks = ax.get_xticks()
integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]
ax.set_xticks(integer_ticks)
ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in integer_ticks])
ax.set_xlabel("$E_{rel}$ Error")
ax.legend(title="Mechanism")

plt.tight_layout()
plt.savefig("error_distribution.png",dpi=300)
plt.show()
