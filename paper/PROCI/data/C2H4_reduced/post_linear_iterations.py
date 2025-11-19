#!python3
import numpy as np
import matplotlib.pyplot as plt
formulations = ["cgc", "cgt"]
label = ["Conservative","Temperature"]
formulation_color  = ["red","blue"]
preconditioners = ["none","jacobi","gauss_seidel"]
preconditioner_labels = ["$P$: none","$P$: Jacobi","$P$: Guass-Seidel"]
preconditioner_style  = ["-","-^","-s"]
for i, formulation in enumerate(formulations):
    for k, preconditioner in enumerate(preconditioners):
        time = np.load(f"{preconditioner}/time_{formulation}.npy")
        iterations = np.load(f"{preconditioner}/linear_iterations_{formulation}.npy")
        plt.plot(time, iterations,  preconditioner_style[k], color = formulation_color[i], label = f"{label[i]} Formulation\n {preconditioner_labels[k]}", markevery=3, mfc = "white")
plt.legend()
plt.xlabel("time [$\mu$ s]")
plt.ylabel("Total number of linear iterations")
plt.xlim([0,15])
plt.show()


