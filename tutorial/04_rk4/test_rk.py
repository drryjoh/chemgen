import numpy as np
import matplotlib.pyplot as plt

#solve dy/dt = a * y
#solution y(t) = y_0 * exp( a * t)

def source(y):
    a = 3
    return a * y

def rk_step(y_n, h):
    k1  = source(y_n)
    k2  = source(y_n +  h * k1/2.0)
    k3  = source(y_n + h * k2/2.0)
    k4  = source(y_n + h * k3)
    return y_n + h/6 * (k1 + 2 * k2 + 2 * k3 + k4)

y0 = 1
max_t = 3
y_n = y0

h1 = 0.1
h2 = 0.01

n1 = int(max_t/h1)
n2 = int(max_t/h2)

sol1 = np.zeros(n1+1)
sol2 = np.zeros(n2+1)
t1 = np.linspace(0,max_t, n1+1)
t2 = np.linspace(0,max_t, n2+1)

sol1[0] = y0
sol2[0] = y0

for i in range(n1):
    sol1[i+1] = rk_step(sol1[i], h1)

for i in range(n2):
    sol2[i+1] = rk_step(sol2[i], h2)

plt.plot(t2, np.exp(3 * t2), '-r', lw = 3)
plt.plot(t1, sol1, '--k', lw = 1)
plt.plot(t2, sol2, '--b', lw = 1)
d = np.loadtxt("out.txt", delimiter=',')
plt.plot(d[:,0], d[:,1], '--g', lw = 1)
plt.show()
