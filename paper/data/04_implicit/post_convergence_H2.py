#!python3
import numpy as np
import matplotlib.pyplot as plt

def L2(data):
    return np.sqrt(np.sum(data**2))

titles = ["SDIRK-4","Rosenbroc","YASS"]
colors = ["red","green", "blue"]
plt.figure(figsize=(8,6))
for k, data_name in enumerate(["sdirk4","rosenbroc","yass"]):
    # Load data
    Ts = np.load(f"{data_name}_Tconvergence.npy")
    Cs = np.load(f"{data_name}_Cconvergence.npy")

    # Compute L2 errors relative to the most refined solution
    print(Ts)
    L2T = np.abs(Ts - Ts[-1])
    CL2 = np.array([L2(C - Cs[-1]) for C in Cs])

    # Print ratios to verify convergence
    print("T convergence ratios:")
    for i in range(len(L2T)-2):
        #print(L2T[i+1])
        print(f"Ref {i+1} -> {i+2}: {L2T[i]/L2T[i+1]}")

    print("\nC convergence ratios:")
    for i in range(len(CL2)-2):
        #print(L2T[i+1])
        print(f"Ref {i+1} -> {i+2}: {CL2[i]/CL2[i+1]}")

    # Prepare plot
    refinement_levels = np.arange(1, len(L2T)+1)
    
    refinement_levels = refinement_levels[:-1]
    L2T = L2T[:-1]
    CL2 =CL2[:-1]

    # Plot L2 errors
    #plt.semilogy(refinement_levels, L2T, 'o-', label='T error')
    plt.semilogy(refinement_levels, CL2, 's-', color=colors[k], label=f'$C$ $L_2$-norm {titles[k]}', mfc="white")

    # Reference lines for T (anchored at L2T[0])
    #plt.semilogy(refinement_levels, L2T[0] * (0.25)**(refinement_levels - 1), '--', label='2nd-order ref (T)')
    #if data_name == "yass":
    #    plt.semilogy(refinement_levels, L2T[0] * (0.5)**(refinement_levels - 1), ':', label='1st-order ref (T)')
    #else:
    #    plt.semilogy(refinement_levels, L2T[0] * (0.0625)**(refinement_levels - 1), ':', label='4th-order ref (T)')
    # Reference lines for C (anchored at CL2[0])
    plt.semilogy(refinement_levels, CL2[0] * (0.25)**(refinement_levels - 1), '--r', color=colors[k], label=f'2nd-order reference')
    if data_name=="yass":
        plt.semilogy(refinement_levels, CL2[0] * (0.5)**(refinement_levels - 1), ':g', color=colors[k], label=f'1st-order reference')
    else:
        plt.semilogy(refinement_levels, CL2[0] * (0.0625)**(refinement_levels - 1), ':g', color=colors[k], label=f'4th-order reference')

    # Formatting
    plt.xlabel("Refinement Level")
    plt.ylabel(f"$C$ $L_2$-norm")
    #plt.title("Convergence of Temperature and Species")
    #plt.xticks(range(1, len(CL2) + 1))  # <- only show 1, 2, 3, 4, ...
    refinement_levels = np.arange(len(CL2))  # 0, 1, 2, 3, ...
    tick_labels = [r"$\Delta t$" if i == 0 else fr"$\Delta t/{2**i}$" for i in refinement_levels]
    plt.xticks(refinement_levels + 1, tick_labels)  # Shift +1 if your levels are labeled 1, 2, 3, 4


    # Only major log grid lines on y-axis
    plt.grid(True, which='major', axis='y')
    plt.minorticks_off()

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3)
plt.tight_layout()
plt.show()
