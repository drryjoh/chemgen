#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)
# Construct file path from directory
L2s = ["errors/L2_FFCM2_model.npy","errors/L2_burke.npy" ]
L2snei = ["errors/L2_nei_FFCM2_model.npy", "errors/L2_nei_burke.npy"]

colors = ["purple", "blue", "orange", "red", "green"]
mechs = ["OConnaire", "burke", "gri30", "sandiego", "FFCM2_model"]
L2s = [f"errors/L2_{mech}.npy" for mech in mechs]
L2snei = [f"errors/L2_nei_{mech}.npy" for mech in mechs]
names = ["Ó Connaire", "Burke", "GRI v3.0", "UCSD", "FFCM2"]

ax = axes[0]
for k, L2 in enumerate(L2s):
    errors = np.load(L2)
    log_errors = np.log10(errors)
    # Plot histogram
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

    # Plot mean line
    mean_val = np.mean(log_errors)
    std_val = np.std(log_errors)
    ax.axvline(mean_val, color=colors[k], linestyle='--', linewidth=1)
xticks = ax.get_xticks()
integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]
ax.set_xticks(integer_ticks)
ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in integer_ticks])

ax = axes[1]
for k, L2 in enumerate(L2snei):
    errors = np.load(L2)
    log_errors = np.log10(errors)
    # Plot histogram
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

    # Plot mean line
    mean_val = np.mean(log_errors)
    std_val = np.std(log_errors)
    ax.axvline(mean_val, color=colors[k], linestyle='--', linewidth=1)
xticks = ax.get_xticks()
integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]
ax.set_xticks(integer_ticks)
ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in integer_ticks])

plt.show()
