#!python3
import numpy as np
import matplotlib.pyplot as plt

#eq 1.3 from ODE II Book
def source_1pt3(y_n):
    s = np.zeros(len(y_n))
    s[0] = -0.04 * y[0] + 1e4 * y[1] * y[2]
    s[1] =  0.04 * y[1] - 1e4 * y[1] * y[2] - 3e7*y[1]**2
    s[2] =  3e7*y[1]**2

def dsource_1pt3_dy(y_n):
    s = np.zeros(len(y_n))
    #s[0] = -0.04 * y[0] + 1e4 * y[1] * y[2]
    s[0][0] = -0.04 
    s[0][1] =  1e4 * y[2] 
    s[0][2] =  1e4 * y[1]

    #s[1] =  0.04 * y[1] - 1e4 * y[1] * y[2] - 3e7*y[1]**2
    s[1] =  0.04 * y[1] - 1e4 * y[1] * y[2] - 3e7*y[1]**2
    s[1] =  0.04 * y[1] - 1e4 * y[1] * y[2] - 3e7*y[1]**2
    s[1] =  0.04 * y[1] - 1e4 * y[1] * y[2] - 3e7*y[1]**2
    
    #s[2] =  3e7*y[1]**2
    s[2] =  3e7*y[1]**2
    s[2] =  3e7*y[1]**2
    s[2] =  3e7*y[1]**2
    

manufactured_x = np.matrix([3.5, 6.6, 8.8, 9.9, 1.1, 0.001, 3.8]).transpose()
manufactured_A = np.matrix(np.zeros([np.shape(manufactured_x)[0],np.shape(manufactured_x)[0]]))
b = manufactured_A * manufactured_x
print("Example 1, A*x")
print(f"A = \n{manufactured_A}")
print(f"y = \n{manufactured_x}")
print(f"b = \n{manufactured_b}")

from scipy.sparse.linalg import gmres

# A: your matrix (dense or sparse)
# b: right-hand side vector

x, info = gmres(A, b)


