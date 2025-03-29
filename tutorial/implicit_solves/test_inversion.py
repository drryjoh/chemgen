#!python3
import numpy as np
import matplotlib.pyplot as plt

#eq 1.3 from ODE II Book
def source_1pt3(y_n):
    y_n = np.asarray(y_n).flatten()  
    s = np.matrix(np.zeros(3)).transpose()
    s[0] = -0.04 * y_n[0] + 1e4 * y_n[1] * y_n[2]
    s[1] =  0.04 * y_n[1] - 1e4 * y_n[1] * y_n[2] - 3e7*y_n[1]**2
    s[2] =  3e7*y_n[1]**2
    return s

def dsource_1pt3_dy(y_n):
    y_n = np.asarray(y_n).flatten()  
    s = np.matrix(np.zeros([3,3]))
    #s[0] = -0.04 * y_n[0] + 1e4 * y_n[1] * y_n[2]
    s[0,0] = -0.04 
    s[0,1] =  1e4 * y_n[2] 
    s[0,2] =  1e4 * y_n[1]

    #s[1] =  0.04 * y_n[1] - 1e4 * y_n[1] * y_n[2] - 3e7*y_n[1]**2
    s[1,0] =  0.0
    s[1,1] =  0.04 - 1e4 * y_n[2] - 2*3e7*y_n[1]
    s[1,2] =  1e4 * y_n[1] 
    
    #s[2] =  3e7*y_n[1]**2
    s[2,0] =  0.0
    s[2,1] =  2*3e7*y_n[1]
    s[2,2] =  0.0
    return s

def derivative_checker(y_n, function, delta):
    #dSdy = dsource_1pt3_dy(y_n)
    #dSdy_check = derivative_checker(y_n, source_1pt3, delta)
    #print(dSdy)
    #print(dSdy_check)
    y_n = np.asarray(y_n).flatten()  
    jacobian = np.zeros([3, 3])
    for i, y_n_i in enumerate(y_n):
        y_n_i_perturbation = np.zeros(3) 
        y_n_i_perturbation[i] = delta
        jacobian[:,i] = (function(y_n + y_n_i_perturbation)-function(y_n))/delta
    return jacobian




#(y^n+1-y^n)/dt - s(y^n+1) = 0
y0 = np.matrix([1,0,0]).transpose()
yn = y0
y_guess_np1 = yn
n_newton = 5
dt = 1e-3
time_final = 0.3
n_time_steps = int(time_final/dt)
time_discretization_jacobian = 1/dt * np.eye(len(yn))

from scipy.sparse.linalg import gmres
ys = [y0]
time = [0]
for t in range(n_time_steps):
    y_guess_np1 = yn + np.matrix([0.0001,0.0001,0.0001]).transpose()
    for k in range(n_newton):
        residual  = (y_guess_np1 - yn)/dt - source_1pt3(y_guess_np1)
        jacobian =  time_discretization_jacobian - dsource_1pt3_dy(y_guess_np1)
        dy, info = gmres(jacobian, -residual)
        dy = np.matrix(dy).transpose()
        y_guess_np1 = y_guess_np1 + dy
    yn = y_guess_np1
    time.append(time[-1]+dt)
    ys.append(yn)
ys = np.array(ys)
time = np.array(time)

plt.plot(time, ys[:,1])

plt.show()







#from scipy.sparse.linalg import gmres

## A: your matrix (dense or sparse)
## b: right-hand side vector

#x, info = gmres(A, b)


