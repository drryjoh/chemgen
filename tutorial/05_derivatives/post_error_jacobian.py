import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data and discard zeros
filename = 'error.txt'
errors = np.loadtxt(filename)
errors = errors[errors > 0]  # Remove zeros

# Set plot style
sns.set(style="whitegrid")

# Plot a histogram on a log scale (log-binned histogram)
log_bins = np.logspace(np.log10(np.min(errors)), np.log10(np.max(errors)), 50)

plt.figure(figsize=(8, 5))
plt.hist(errors, bins=log_bins, edgecolor='black')
plt.xscale('log')

plt.title('Log-Scale Error Frequency Histogram')
plt.xlabel('Error (log scale)')
plt.ylabel('Frequency')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()
