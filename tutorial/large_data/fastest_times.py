#!python3
import numpy as np
import matplotlib.pyplot as plt

# Load the eigenvalues
data = []
fastest = []
n_fast = 0

# Load eigenvalues
eigenvalues = np.load("eigenvalues.npy")
f2 = open("fastest.csv", "w")

with open("consolidated.csv", "r") as f:
    for value in eigenvalues:
        try:
            line = next(f)
        except StopIteration:
            print("Warning: Reached end of file before reading all eigenvalues.")
            break
        if np.log10(value) < -7.7:
            f2.write(line)  # strip newline characters
            fastest.append(np.log10(value))
            n_fast += 1
print(n_fast)

#Remove the placeholder values (1.0)
filtered = np.log10(eigenvalues[eigenvalues != 1.0])

# Plot histogram
plt.figure(figsize=(10, 6))

# Compute shared bins
bins = np.histogram_bin_edges(np.concatenate((filtered, fastest)), bins=100)

# Plot both histograms using the same bins
plt.hist(filtered, bins=bins, edgecolor='k', label='All')
plt.hist(fastest, bins=bins, edgecolor='k', color='red', alpha=0.5, label='Fastest')

plt.legend()
plt.title("Frequency Distribution of Eigenvalues (T > 1000)")
plt.xlabel("Eigenvalue")
plt.ylabel("Frequency")
plt.grid(True)
plt.tight_layout()
plt.show()

