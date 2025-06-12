#!python3
import numpy as np
import matplotlib.pyplot as plt
# Load ChemGen output
plt.rcParams["lines.linewidth"] = 2.5
d4 = np.loadtxt("chem_out_4.txt")
d7 = np.loadtxt("chem_out_7.txt")
d12 = np.loadtxt("chem_out_12.txt")
ct_T = np.load("c_T.npy")
ct_time = np.load("c_time.npy")

# Set up side-by-side plots

# ---- Plot Temperature Evolution (Left) ----
plt.plot(ct_time*1000.0, ct_T, '-r', label="Cantera")
plt.plot(d4[:, 0]*1000.0, d4[:, 1], '--g', label="ChemGen, $n_p=4$")
plt.plot(d7[:, 0]*1000.0, d7[:, 1], '--k', label="ChemGen, $n_p=7$")
plt.plot(d12[:, 0]*1000.0, d12[:, 1], '--b', label="ChemGen, $n_p=12$")

ax = plt.gca()
ax.set_xlim([0.05, 0.2])
ax.set_ylim([2600,2900])
ax.set_xlabel("Time ($\mu$s)")
ax.set_ylabel("Temperature (K)")
ax.set_title("Temperature Evolution")
plt.legend()
plt.savefig("temperature_np.png",dpi=300)
plt.show()
