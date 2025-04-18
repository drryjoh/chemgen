#!python3
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import os

# Load master data
original_data = np.loadtxt("consolidated.csv", delimiter=',')
T = original_data[:, 0]
n_points = len(T)

# Function to compute eigenvalue or 1.0
def process_matrix(i):
    if T[i] <= 1000:
        return 1.0
    A = np.loadtxt(f"jacobians/jacobian_{i}.txt", delimiter=',')
    eigvals = np.linalg.eigvals(A)
    return 1.0/np.max(np.abs(np.real(eigvals)))  # or any other logic to pick eigenvalue

# Parallel loop
if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        eigen_values = list(executor.map(process_matrix, range(n_points)))

    # Convert to NumPy array
    eigen_values = np.array(eigen_values)
    np.save("eigenvalues.npy", eigen_values)
