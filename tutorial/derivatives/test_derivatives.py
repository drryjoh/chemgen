import numpy as np
import matplotlib.pyplot as plt

def function(x):
    return x**3.4

def dfunction_dx(x):
    return 3.4*x**2.4

def finite_difference(f,x, dx):
    return (f(x+dx)-f(x-dx))/(2*dx)

starting_dx = 3
x = 5
n_refine = 5
dfdx = []
dxs = []

for i in range(n_refine):
    dx = starting_dx * (1/2)**i
    dxs.append(dx)
    dfdx.append(finite_difference(function, x, dx))
dfdx_analytical = dfunction_dx(x)
print(dfdx)
np.array(dxs)
dfdx_er = np.array(np.array(dfdx)-dfdx_analytical)
plt.loglog(dxs, dfdx_er,'-r')
plt.loglog(dxs,  np.array([dfdx_er[0] * (1/2)**(2*i) for i in range(n_refine)]),'--k')
plt.xlabel("dx")
plt.ylabel("Error wrt Analytical Derivative")
plt.show()


