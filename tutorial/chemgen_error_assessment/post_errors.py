import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cantera as ct

# Load the data from your CSV file
data = np.loadtxt("l2_norm_results.csv", delimiter=",", skiprows=1)
states = np.loadtxt("states.csv", delimiter=",", skiprows=0)

# Separate the data into temperatures and errors
temperatures = data[:, 0]  # First column
errors = data[:, 1]        # Second column
bad_states = []
for k, error in enumerate(errors):
    if error>10.0:
        bad_states.append(states[k,:])
        np.save(f"bad_state_{k}.npy", states[k,:])
for bad_state in bad_states:
    gas = ct.Solution("FFCM2_model.yaml")
    temperature = bad_state[0]
    pressure = ct.gas_constant * temperature * np.sum(bad_state[1:])
    X = bad_state[1:]/np.sum(bad_state[1:])
    gas.TPX =  temperature, pressure, X
    print(pressure)
    print(X)
# Create a boxplot
fig, ax = plt.subplots(figsize=(8, 6))
boxprops = dict(patch_artist=True)  # Enable patching for the box
boxplot = ax.boxplot(errors, vert=True, patch_artist=True, showmeans=True, flierprops={'marker': 'o', 'color': 'black', 'alpha': 0},
                     meanprops=dict(marker='D', markeredgecolor='k', markerfacecolor='red'),
                     medianprops=dict(color='k', linewidth=1, linestyle = '--'))


# Log-scale the y-axis
ax.set_yscale('log')

# Normalize temperatures for color mapping
norm = Normalize(vmin=np.min(temperatures), vmax=np.max(temperatures))
cmap = plt.cm.viridis  # Choose a colormap

# Get outlier points (fliers)
for flier in boxplot['fliers']:
    flier.set_marker('o')  # Circle marker for outliers
    flier.set_alpha(0)     # Make the default flier invisible, handled below

# Manually plot colored outliers
outlier_indices = [i for i, e in enumerate(errors) if e < boxplot['whiskers'][0].get_ydata()[0] or e > boxplot['whiskers'][1].get_ydata()[1]]
outlier_errors = errors[outlier_indices]
outlier_temperatures = temperatures[outlier_indices]

scatter = ax.scatter(
    np.full_like(outlier_errors, 1),  # Plot all at x=1
    outlier_errors,
    c=outlier_temperatures,
    cmap=cmap,
    norm=norm,
    edgecolor=None,
    s=20,  # Adjust size of the circles
)

# Add a colorbar
cbar = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax)
cbar.set_label("Temperature", fontsize=14)

# Customize the plot
ax.set_title("Box-and-Whisker Plot with Colored Outliers (Log Scale)", fontsize=16)
ax.set_ylabel("Relative Error (log scale)", fontsize=14)
ax.set_xticks([1])
ax.set_xticklabels(["Errors"])
ax.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

# Show the plot
plt.show()
