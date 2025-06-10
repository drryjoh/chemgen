#!python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cantera as ct
import seaborn as sns

# Load the data from your CSV file
fig, ax = plt.subplots(figsize=(8, 6))
colors = ["purple", "blue", "orange", "red", "green"]   # Match number of mechs
mechs = ["OConnaire","burke","gri30","FFCM2_model","sandiego"]
names = ["FFCM2 "]
for i, mech in enumerate(mechs):
    data = np.loadtxt(f"l2_{mech}.csv", delimiter=",", skiprows=1)
    states = np.loadtxt(f"{mech}.csv", delimiter=",", skiprows=0)
    # Separate the data into temperatures and errors
    temperatures = data[:, 0]  # First column
    errors = data[:, 1]        # Second column
    sns.histplot(np.log10(errors), ax=ax, kde=True, bins = 200, color = colors[i], alpha = 0.4, edgecolor =  'black', label = mech)

ax.legend(title="Mechanism")
xticks = ax.get_xticks()
integer_ticks = xticks[np.isclose(xticks, np.round(xticks))]  # Only integers
ax.set_xticks(integer_ticks)

ax.set_title("Error Distribution", fontsize=16)
ax.set_xlabel("$\epsilon = \\frac{1}{n_s} \sum_{i=1}^{n_s} |(s(y)_{i,cg}-s(y)_{i,ct})/s(y)_{i,ct}|$", fontsize=10)
ax.set_ylabel("Error Distribution", fontsize=14)

# Replace xtick labels with 10^number format
ax.set_xticklabels([f"$10^{{{int(tick)}}}$" for tick in integer_ticks])
plt.tight_layout()
plt.savefig("hist.png", dpi=300)

# Show the plot
plt.show()
