#!python3
import chemgen as cg
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def random_yc(ns):
    species_array =  np.random.uniform(0, 1, ns)
    species_array /= species_array.sum()

    return (1000 + 1500 * np.random.random(), 10132.5 + 101325.0 * 9.9 * np.random.random(), species_array)

def discrete_jacobian(C,T,dc):
    n_species  = len(C)
    J = np.zeros([ns, ns])
    for i in range(n_species):
        C[i] = C[i] + dc
        S1 = np.array(cg.source(C, T))
        C[i] = C[i] - dc
        C[i] = C[i] - dc
        S2 = np.array(cg.source(C, T))
        C[i] = C[i] + dc
        J[:,i] = (S1-S2) / (2 * dc)
    return J
def Frobenius(J, ns):
    L2 = 0
    for i in range(ns):
        for j in range(ns):
            L2 += J[i,j]**2
    return np.sqrt(L2)
    

def check_J(Jcg, Jfd, ns):
    minJ = np.min(Jfd[np.abs(Jfd)>0])/1000
    L2_J = Frobenius(Jcg, ns)
    for i in range(ns):
        for j in range(ns):
            L2 = 0
            if np.abs(Jcg[i,j]) > L2_J/(10**20):
                L2 = np.abs((Jcg[i,j] - Jfd[i,j])/(Jfd[i,j]))
            if L2 > 1e-3:
                print(f"{Jcg[i,j]} {Jfd[i,j]} {(Jcg[i,j] - Jfd[i,j])/(Jcg[i,j]) }")

def L2_J(Jcg, Jfd, ns):
    L2 = 0
    L2_J = Frobenius(Jcg, ns)
    for i in range(ns):
        for j in range(ns):
            if np.abs(Jcg[i,j]) > L2_J/(10**20):
                L2 += ((Jcg[i,j] - Jfd[i,j]))**2
    return np.sqrt(L2)/L2_J

def L2_nei_J(Jcg, Jfd, ns):
    L2 = 0
    L2_J = Frobenius(Jcg, ns)
    number_of_elements = 0
    for i in range(ns):
        for j in range(ns):
            if np.abs(Jcg[i,j]) > L2_J/(10**10):
                L2 += ((Jcg[i,j] - Jfd[i,j])/Jcg[i,j])**2
                number_of_elements+=1
    print(number_of_elements/ns**2)
    return np.sqrt(L2)

mech = "sandiego"
gas = ct.Solution(f"{mech}.yaml")

ns = gas.n_species

#create random chemical state:
n_random = 10000
yc = []

gas.TPX = random_yc(ns)
C = gas.concentrations
T = gas.T
yc.append(np.concatenate(([T], C)))

L2s = []
yci = yc[0]
T = yci[0]
C = yci[1:]
dc = np.min(C[C > 0])/1.5
Jcg = np.array(cg.source_jacobian(C,T))
for r in range(4):
    Jfd = discrete_jacobian(C, T, dc/(2**r))
    L2s.append(L2_J(Jcg, Jfd, ns))
refinement  = np.array(range(4))
L2s = np.array(L2s)

plt.figure()
plt.semilogy(refinement, L2s / L2s[0], 'ok', label='Observed error')
plt.semilogy(refinement, [1 / (4 ** r) for r in refinement], '--r', label='Second-order reference')

# Custom x-axis tick labels
tick_labels = [f"[$\delta c/{2**r}$, $\delta T/{2**r}$]" if r > 0 else "[$\delta c$, $\delta T$]" for r in refinement]
plt.xticks(refinement, tick_labels)

# Label the plot
plt.xlabel('Refinement Level (r)')
plt.ylabel('Normalized L2 Error')
plt.legend()
plt.grid(True, which="both", ls="--", lw=0.5)
plt.savefig("ooa.png",dpi=300)
plt.show()