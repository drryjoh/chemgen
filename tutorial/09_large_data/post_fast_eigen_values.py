#!python3
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import os

# Load master data
original_data = np.loadtxt("fastest.csv", delimiter=',')
T = original_data[:, 0]
n_points = len(T)
n_reactions = 34

# Function to compute eigenvalue or 1.0
def process_matrix(i):
    if i%100:
        print(i)
    dat = []
    for j in range(n_reactions):
        A = np.loadtxt(f"data_fast/R_{i}_{j}.txt", delimiter=',')
        eigvals = np.linalg.eigvals(A)
        dat.append(1.0/np.max(np.abs(np.real(eigvals))))  # or any other logic to pick eigenvalue
    return dat

# Parallel loop
if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        eigen_values = list(executor.map(process_matrix, range(n_points)))

    # Convert to NumPy array
    eigen_values = np.array(eigen_values)
    np.save("eigenvalues_fast.npy", eigen_values)
