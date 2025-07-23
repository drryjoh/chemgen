#!python3
import chemgen as cg
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns

def eps():
    return 100.*np.finfo(float).eps

def random_yc(ns):
    species_array =  np.random.uniform(eps(), 1. - eps(), ns)
    species_array /= species_array.sum()

    return (1000 + 1500 * np.random.random(), 10132.5 + 101325.0 * 9.9 * np.random.random(), species_array)

def discrete_jacobian(yci, dc, dT, temperature_equation, ignore_temp_dependence):
    J = np.zeros([nv, nv])

    if ignore_temp_dependence:
        T = cg.temperature_from_internal_energy(yci[1:], yci[0])
        for i in range(nv):
            if i == 0:
                dy = dT
            else:
                dy = dc[i-1]
            yci[i] = yci[i] + dy
            S1 = np.array(cg.source(yci[1:], T))
            yci[i] = yci[i] - 2.*dy
            S2 = np.array(cg.source(yci[1:], T))
            yci[i] = yci[i] + dy
            J[:,i] = (S1-S2) / (2 * dy)
    else:
        for i in range(nv):
            if i == 0:
                dy = dT
            else:
                dy = dc[i-1]
            yci[i] = yci[i] + dy
            T1 = cg.temperature(yci[1:], yci[0])
            S1 = np.array(cg.source(yci[1:], T1))
            yci[i] = yci[i] - 2.*dy
            T2 = cg.temperature(yci[1:], yci[0])
            S2 = np.array(cg.source(yci[1:], T2))
            yci[i] = yci[i] + dy
            J[:,i] = (S1-S2) / (2 * dy)

    return J

def Frobenius(J, nv):
    L2 = 0
    for i in range(nv):
        for j in range(nv):
            L2 += J[i,j]**2
    return np.sqrt(L2)


# def check_J(Jcg, Jfd, ns):
#     minJ = np.min(Jfd[np.abs(Jfd)>0])/1000
#     L2_J = Frobenius(Jcg, ns)
#     for i in range(ns):
#         for j in range(ns):
#             L2 = 0
#             if np.abs(Jcg[i,j]) > L2_J/(10**20):
#                 L2 = np.abs((Jcg[i,j] - Jfd[i,j])/(Jfd[i,j]))
#             if L2 > 1e-3:
#                 print(f"{Jcg[i,j]} {Jfd[i,j]} {(Jcg[i,j] - Jfd[i,j])/(Jcg[i,j]) }")

def L2_J(Jcg, Jfd, nv):
    L2 = 0
    L2_J = Frobenius(Jcg, nv)
    Jdiff = (Jcg - Jfd)/(Jcg + eps())
    # breakpoint()
    for i in range(nv):
        for j in range(nv):
            # if i != 0 and j != 0: continue # uncomment to check only derivatives associated with energy variable
            # if i == 0 or j == 0: continue # uncomment to check only derivatives associated with concentrations
            # if not (i == 0 and j == 0): continue
            # if not (i >= 0 and j == 0): continue
            # if not (j >= 0 and i == 0): continue
            # print("({i}, {j})".format(i=i, j=j))
            if np.abs(Jcg[i,j]) > L2_J/(10**20):
                L2 += Jdiff[i,j]**2
    return np.sqrt(L2)

# def L2_nei_J(Jcg, Jfd, nv):
#     L2 = 0
#     L2_J = Frobenius(Jcg, nv)
#     number_of_elements = 0
#     for i in range(nv):
#         for j in range(nv):
#             if np.abs(Jcg[i,j]) > L2_J/(10**10):
#                 L2 += ((Jcg[i,j] - Jfd[i,j])/Jcg[i,j])**2
#                 number_of_elements+=1
#     print(number_of_elements/nv**2)
#     return np.sqrt(L2)

mech = "ffcm2_h2"
gas = ct.Solution(f"{mech}.yaml")

ns = gas.n_species
nv = ns + 1

#create random chemical state:
yc = []

TPX = random_yc(ns)
# TPX = (1455.3201340478888, 852156.9373644177, np.array([0.00417861, 0.11311578, 0.85521151, 0.0274941 ]))
# TPX = (1947.922280896753, 957477.6343142525, np.array([0.09818607, 0.22232688, 0.0175851 , 0.0890133 , 0.00461794,
#        0.19607955, 0.06379871, 0.12727434, 0.18111811]))
gas.TPX = TPX
C = gas.concentrations
T = gas.T
yc.append(np.concatenate(([T], C)))

L2s = []
yci = yc[0]
# T = yci[0]
# C = yci[1:]

if cg.temperature_equation():
    yci[0] = T
else:
    yci[0] = cg.internal_energy_volume_specific(C, T)

# dc = np.min(C[C > 0])/1.5
dc = C/1.e1
dT = yci[0]/1.e3
Jcg = np.array(cg.source_jacobian(C,T))
# breakpoint()
for r in range(4):
    Jfd = discrete_jacobian(yci, dc/(2**r), dT/(2**r), cg.temperature_equation(), cg.ignore_temp_dependence())
    L2s.append(L2_J(Jcg, Jfd, nv))
refinement  = np.array(range(4))
L2s = np.array(L2s)

# breakpoint()
print("TPX = {}".format(TPX))
print("dc = {}".format(dc))
print("dT = {}".format(dT))
print("L2s = {}".format(L2s))

plt.figure()
plt.semilogy(refinement, [1.0 / (4 ** r) for r in refinement], '--r', label='Second-order reference', lw=3)
plt.semilogy(refinement, L2s / L2s[0], ':ok', label='Observed $L_2$')


# Custom x-axis tick labels
tick_labels = [f"[$\\delta c/{2**r}$, $\\delta T/{2**r}$]" if r > 0 else "[$\\delta c$, $\\delta T$]" for r in refinement]
plt.xticks(refinement, tick_labels, fontsize=12)

# Label the plot
plt.xlabel('Refinement Level',fontsize=16)
plt.ylabel('$L_2$',fontsize=16)
plt.legend(fontsize=16)
plt.grid(True, which="both", ls="--", lw=0.5)
plt.savefig("ooa.png",dpi=300)
plt.show()