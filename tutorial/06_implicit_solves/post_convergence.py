#!python3
import numpy as np
import matplotlib.pyplot as plt

def L2(data):
    return np.sqrt(np.sum(data**2))

# Load data
Ts = np.load("Tconvergence.npy")
Cs = np.load("Cconvergence.npy")

# Compute L2 errors relative to the most refined solution
L2T = np.abs(Ts - Ts[-1])
CL2 = np.array([L2(C - Cs[-1]) for C in Cs])

# Print ratios to verify convergence
print("T convergence ratios:")
for i in range(len(L2T)-1):
    print(f"Ref {i+1} -> {i+2}: {L2T[i]/L2T[i+1]:.2f}")

print("\nC convergence ratios:")
for i in range(len(CL2)-1):
    print(f"Ref {i+1} -> {i+2}: {CL2[i]/CL2[i+1]:.2f}")

# Prepare plot
refinement_levels = np.arange(1, len(L2T)+1)
plt.figure(figsize=(8,6))
refinement_levels = refinement_levels[:-1]
L2T = L2T[:-1]
CL2 =CL2[:-1]

# Plot L2 errors
plt.semilogy(refinement_levels, L2T, 'o-', label='T error')
plt.semilogy(refinement_levels, CL2, 's-', label='C error')

# Reference lines for T (anchored at L2T[0])
plt.semilogy(refinement_levels, L2T[0] * (0.25)**(refinement_levels - 1), '--', label='2nd-order ref (T)')
plt.semilogy(refinement_levels, L2T[0] * (0.0625)**(refinement_levels - 1), ':', label='4th-order ref (T)')

# Reference lines for C (anchored at CL2[0])
plt.semilogy(refinement_levels, CL2[0] * (0.25)**(refinement_levels - 1), '--', label='2nd-order ref (C)')
plt.semilogy(refinement_levels, CL2[0] * (0.0625)**(refinement_levels - 1), ':', label='4th-order ref (C)')

# Formatting
plt.xlabel("Refinement Level")
plt.ylabel("L2 Error")
plt.title("Convergence of Temperature and Species")
plt.legend()

# Only major log grid lines on y-axis
plt.grid(True, which='major', axis='y')
plt.minorticks_off()

plt.tight_layout()
plt.show()
