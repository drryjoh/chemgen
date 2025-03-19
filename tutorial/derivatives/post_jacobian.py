import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

def plot_heatmap(matrix):
    plt.figure(figsize=(8, 6))

    # Set the normalization to SymLogNorm to handle negative and small values
    norm = mcolors.SymLogNorm(linthresh=1e3, linscale=1, vmin=np.min(matrix), vmax=np.max(matrix))

    sns.heatmap(matrix, annot=True, cmap="viridis", fmt="3.2e", linewidths=0.5, 
                cbar=True, annot_kws={"size": 4}, norm=norm)

    plt.title("Jacobian at Intermediate State")
    plt.xlabel("Columns")
    plt.ylabel("Rows")
    plt.show()

# Load the matrix from file
matrix = np.loadtxt("jacobian_out.txt", delimiter=',')
eigenvalues = np.linalg.eigvals(matrix)
print("Eigenvalues of the matrix:")
print(eigenvalues)

plot_heatmap(matrix)
