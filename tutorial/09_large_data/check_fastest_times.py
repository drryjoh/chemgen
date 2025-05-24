#!python3
import numpy as np
import matplotlib.pyplot as plt

# Load the eigenvalues
data = []
fastest = []
n_fast = 0
n_reactions = 34

# Load eigenvalues
eigenvalues = np.load("eigenvalues.npy")
filtered = eigenvalues[eigenvalues != 1.0]
R = np.load("eigenvalues_fast.npy")

for j in range(n_reactions):
    fastest = R[:,j]
    # Plot histogram
    plt.figure(figsize=(10, 6))

    # Compute shared bins
    # Get min and max of combined data
    combined = np.concatenate((fastest, filtered))
    log_min = np.log10(np.min(combined))
    log_max = np.log10(np.max(combined))

    # Create log-spaced bins
    bins = np.logspace(log_min, log_max, num=100)

    # Plot both histograms using the same bins
    plt.hist(filtered, bins=bins, edgecolor='k', label='All')
    plt.hist(fastest, bins=bins, edgecolor='k', color='red', alpha=0.5, label='Fastest')
    plt.xscale('log')

    plt.legend()
    plt.title(f"Frequency Distribution of Eigenvalues R = {j}")
    plt.xlabel("Eigenvalue")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"images/R{j}.png",dpi=300)


