#!python3
import numpy as np
import matplotlib.pyplot as plt

# Load the eigenvalues
eigenvalues = np.load("eigenvalues.npy")

# Remove the placeholder values (1.0)
filtered = eigenvalues[eigenvalues != 1.0]

# Plot histogram
plt.figure(figsize=(10, 6))
log_min = np.log10(min(filtered))
log_max = np.log10(max(filtered))
bins = np.logspace(log_min, log_max, num=100)
plt.hist(filtered, bins=bins, edgecolor='k')
plt.title("Frequency Distribution of Eigenvalues (T > 1000)")

plt.xlabel("Eigenvalue")
plt.xscale('log')
plt.ylabel("Frequency")
plt.grid(True)
plt.tight_layout()
plt.savefig("histogram_eigenvalues.png",dpi=300)

