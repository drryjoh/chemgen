import numpy as np
import sys
import matplotlib.pyplot as plt

mechanisms  = ["chemgen_ffcm2", "chemgen_h2"]
try:
    ptsa = np.load(f"data/{mechanisms[0]}_pts.npy")
except FileNotFoundError:
    print(f"Error: The file 'data/{mechanisms[0]}_pts.npy' was not found.")
    print("Please run 'grab_sherlock_data.py' to generate the required file.")
    sys.exit(1)

ptsb = np.load(f"data/{mechanisms[1]}_pts.npy")
kfa_threaded = np.load(f"data/{mechanisms[0]}_forward_reactions_nthreads.npy")
kfb_threaded = np.load(f"data/{mechanisms[1]}_forward_reactions_nthreads.npy")
kfa_serial = np.load(f"data/{mechanisms[0]}_forward_reactions_serial.npy")
kfb_serial = np.load(f"data/{mechanisms[1]}_forward_reactions_serial.npy")

proga_threaded = np.load(f"data/{mechanisms[0]}_point_progress_rates_nthreads.npy")
progb_threaded = np.load(f"data/{mechanisms[1]}_point_progress_rates_nthreads.npy")
proga_serial = np.load(f"data/{mechanisms[0]}_point_progress_rates_serial.npy")
progb_serial = np.load(f"data/{mechanisms[1]}_point_progress_rates_serial.npy")

sourcea_threaded = np.load(f"data/{mechanisms[0]}_source_nthreads.npy")
sourceb_threaded = np.load(f"data/{mechanisms[1]}_source_nthreads.npy")
sourcea_serial = np.load(f"data/{mechanisms[0]}_source_serial.npy")
sourceb_serial = np.load(f"data/{mechanisms[1]}_source_serial.npy")

plt.semilogx(ptsa, kfa_serial/kfa_threaded, '-ob', mfc="white")
plt.semilogx(ptsb, kfb_serial/kfb_threaded, '-or', mfc="white")

plt.semilogx(ptsa, proga_serial/proga_threaded, '-^b', mfc="white")
plt.semilogx(ptsb, progb_serial/progb_threaded, '-^r', mfc="white")

plt.semilogx(ptsa, sourcea_serial/sourcea_threaded, '-xb', mfc="white")
plt.semilogx(ptsb, sourceb_serial/sourceb_threaded, '-xr', mfc="white")

# Proxy artists for the legend
lines = [
    # Column 1: Line descriptions
    plt.Line2D([], [], color='blue', linestyle='-', label="Full FFCM2 Mechanism,\n$n_s=96$ and $n_r=1073$"),
    plt.Line2D([], [], color='red', linestyle='-', label="FFCM2 H2 Mechanism,\n$n_s=12$ and $n_r=32$"),
    plt.Line2D([], [], color='none', linestyle='None', label=""),  # Spacer
    plt.Line2D([], [], color='black', marker='o', mfc="white", linestyle='None', label="Forward rates, $k_{f,j}$"),
    plt.Line2D([], [], color='black', marker='^', mfc="white", linestyle='None', label="Progress rates, $q_{j}$"),
    plt.Line2D([], [], color='black', marker='x', mfc="white", linestyle='None', label="Source term, $\\omega_{i}$"),
]

# Create a two-column legend
plt.legend(handles=lines, loc="best", ncol=2, handletextpad=1.5, columnspacing=2, fontsize=8)
plt.ylim([0,15])
plt.yticks([0,5,10,15])
plt.xlabel("Number of points")
plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")

plt.figure()
besta = np.load(f"data/{mechanisms[0]}_best_time.npy")
seriala = np.load(f"data/{mechanisms[0]}_serial_time.npy")
bestb = np.load(f"data/{mechanisms[1]}_best_time.npy")
serialb = np.load(f"data/{mechanisms[1]}_serial_time.npy")

plt.semilogx(ptsa, seriala/besta,'-ob', label="Full FFCM2 Mechanism,\n $n_s=96$ and $n_r=1073$")
plt.semilogx(ptsb, serialb/bestb,'-or', label="FFCM2 H2 Mechanism,\n $n_s=12$ and $n_r=32$")
plt.legend()
plt.xlabel("Number of points")
plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")
# Show the plot (for demonstration purposes)
plt.show()


