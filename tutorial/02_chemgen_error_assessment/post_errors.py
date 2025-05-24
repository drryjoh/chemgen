#!python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cantera as ct
import seaborn as sns

# Load the data from your CSV file
data = np.loadtxt("l2_norm_results.csv", delimiter=",", skiprows=1)
states = np.loadtxt("states.csv", delimiter=",", skiprows=0)

# Separate the data into temperatures and errors
temperatures = data[:, 0]  # First column
errors = data[:, 1]        # Second column
bad_states = []
for k, error in enumerate(errors):
    if error>.8:
        bad_states.append(states[k,:])
        np.save(f"bad_state_{k}.npy", states[k,:])
for bad_state in bad_states:
    gas = ct.Solution("FFCM2_model.yaml")
    temperature = bad_state[0]
    pressure = ct.gas_constant * temperature * np.sum(bad_state[1:])
    X = bad_state[1:]/np.sum(bad_state[1:])
    gas.TPX =  temperature, pressure, X
    print(pressure)
    for Xi, species in zip(X, gas.species_names):
        if Xi>0.0:
            print(f"{species}: {Xi}")
# Create a boxplot
fig, ax = plt.subplots(figsize=(8, 6))
boxprops = dict(patch_artist=True)  # Enable patching for the box
boxplot = ax.boxplot(errors, vert=True, patch_artist=True, showmeans=True, flierprops={'marker': 'o', 'color': 'black', 'alpha': 0},
                     meanprops=dict(marker='D', markeredgecolor='k', markerfacecolor='red'),
                     medianprops=dict(color='k', linewidth=1, linestyle = '--'))
for i, flier in enumerate(boxplot['fliers']):
    outliers = flier.get_ydata()  # Get the y-coordinates of the outliers
    x_pos = np.full(outliers.shape, i + 1)  # X-position for the box (1-based index)
    ax.plot(x_pos, outliers, 'o', color='k', mfc='white', markersize = 4, label="Outliers" if i == 0 else "")  # Plot outliers as points


# Log-scale the y-axis
ax.set_yscale('log')


# Customize the plot

ax.set_ylabel("Relative Error L2-norm\n $L_2 = \sqrt{\sum_{i=1}^{n_s}\\frac{1}{n_s}((s(y)_{i,cg}-s(y)_{i,ct})/s(y)_{i,ct})^2}$", fontsize=10)
ax.set_xticks([1])
ax.set_xticklabels(["Box-whisker plot with outliers"])
ax.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig("bw.png", dpi=300)

# Distribution plot (KDE or histogram)
fig, ax = plt.subplots(figsize=(8, 6))
sns.histplot(np.log10(errors), kde=True, ax=ax, bins=1000, color='k', alpha=0.6, edgecolor='black')
ax.set_title("Error Distribution", fontsize=16)
ax.set_xlabel("Relative Error L2-norm\n $L_2 = \sqrt{\sum_{i=1}^{n_s}\\frac{1}{n_s}((s(y)_{i,cg}-s(y)_{i,ct})/s(y)_{i,ct})^2}$", fontsize=10)
ax.set_ylabel("Frequency", fontsize=14)
xticks = ax.get_xticks()
integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]  # Only integers
ax.set_xticks(integer_ticks)
plt.tight_layout()


# Replace xtick labels with 10^number format
ax.set_xticklabels([f"$10^{{{int(tick)}}}$" for tick in integer_ticks])
plt.savefig("hist.png", dpi=300)

# Show the plot
plt.show()
