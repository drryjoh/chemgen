import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

def plot_heatmap(matrix, file, i):
    plt.figure(figsize=(8, 6))

    # Set the normalization to SymLogNorm to handle negative and small values
    norm = mcolors.SymLogNorm(linthresh=1e3, linscale=1, vmin=np.min(matrix), vmax=np.max(matrix))

    sns.heatmap(matrix, annot=True, cmap="viridis", fmt="3.2e", linewidths=0.5, 
                cbar=True, annot_kws={"size": 4}, norm=norm)

    plt.title("Jacobian at Intermediate State")
    plt.xlabel("Columns")
    plt.ylabel("Rows")
    plt.savefig(f"images/{file}_{i}.png")

# Load the matrix from file
file = "h2o2"
n_time_steps = 300
stiffness = np.zeros(n_time_steps)
condition = np.zeros(n_time_steps)
for i in range(n_time_steps):
    matrix = np.loadtxt(f"data/jacobian_out_{file}_{i}.txt", delimiter=',')
    eigenvalues = np.linalg.eigvals(matrix)
    real_eigenvalues = eigenvalues.real
    abs_nonzero_real_eigenvalues = np.abs(real_eigenvalues[real_eigenvalues != 0])
    stiffness[i] =  np.max(abs_nonzero_real_eigenvalues)/np.min(abs_nonzero_real_eigenvalues)
# Step 2: Compute Condition Number
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    condition_number = max(singular_values) / min(singular_values)
    condition[i] = condition_number
    plot_heatmap(matrix, file, i)
    print(f"Calculated Eigen Values for Step {i}")
np.save(f"stiffness_{file}.npy",stiffness)
np.save(f"condition_{file}.npy",condition)
