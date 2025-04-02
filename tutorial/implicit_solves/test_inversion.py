#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import time

import numpy as np

def back_substitution(R, b):
    """
    Solves Rx = b where R is upper triangular.
    """
    n = len(b)
    x = np.zeros_like(b)
    for i in reversed(range(n)):
        x[i] = (b[i] - np.dot(R[i, i+1:], x[i+1:])) / R[i, i]
    return x

def least_squares_qr(A, b):
    """
    Solves min ||Ax - b|| using QR decomposition (full-rank).
    Equivalent to np.linalg.lstsq(A, b)[0]
    """
    Q, R = np.linalg.qr(A)
    Qt_b = Q.T @ b
    x = back_substitution(R, Qt_b)
    return x

def norm2(x):
    return np.sqrt(np.sum(x**2))

def gmres_custom(A, b, x0=None, tol=1e-10, max_iter=100):
    """
    Minimal GMRES solver for Ax = b.

    Parameters:
        A        : function or matrix (2D ndarray or callable that applies A to a vector)
        b        : right-hand side vector
        x0       : initial guess (default is zero)
        tol      : residual tolerance
        max_iter : max number of iterations

    Returns:
        x        : approximate solution
        info     : 0 if converged, 1 otherwise
    """
    n = len(b)
    if x0 is None:
        x0 = np.zeros_like(b)

    def matvec(v):
        return A @ v if callable(getattr(A, '__matmul__', None)) else A(v)

    r0 = b - matvec(x0)
    beta = norm2(r0)
    if beta < tol:
        return x0, 0

    V = np.zeros((n, max_iter + 1))
    H = np.zeros((max_iter + 1, max_iter))
    V[:, 0] = r0 / beta
    g = np.zeros(max_iter + 1)
    g[0] = beta

    for j in range(max_iter):
        w = matvec(V[:, j])

        # Modified Gram-Schmidt
        for i in range(j + 1):
            H[i, j] = np.dot(V[:, i], w)
            w -= H[i, j] * V[:, i]

        H[j + 1, j] = norm2(w)
        if H[j + 1, j] != 0 and j + 1 < n:
            V[:, j + 1] = w / H[j + 1, j]

        # Solve the least squares problem min ||g - H y||
        y = least_squares_qr(H[:j + 2, :j + 1], g[:j + 2])
        x = x0 + V[:, :j + 1] @ y
        residual = norm2(matvec(x) - b)

        if residual < tol:
            return x, 0

    return x, 1  # Did not converge


def timeit_solver(name, func, *args, **kwargs):
    print(f"Running {name}...")
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    print(f"{name} completed in {elapsed:.4f} seconds\n")
    return result, elapsed


def source_1pt3(y_n):
    y_n = np.asarray(y_n).flatten()
    s = np.zeros(3)
    s[0] = -0.04 * y_n[0] + 1e4 * y_n[1] * y_n[2]
    s[1] =  0.04 * y_n[0] - 1e4 * y_n[1] * y_n[2] - 3e7 * y_n[1]**2
    s[2] =  3e7 * y_n[1]**2
    return s

def dsource_1pt3_dy(y_n):
    y_n = np.asarray(y_n).flatten()
    s = np.zeros((3,3))
    #s[0] = -0.04 * y_n[0] + 1e4 * y_n[1] * y_n[2]
    s[0,0] = -0.04 
    s[0,1] =  1e4 * y_n[2] 
    s[0,2] =  1e4 * y_n[1]
    #s[1] =  0.04 * y_n[0] - 1e4 * y_n[1] * y_n[2] - 3e7 * y_n[1]**2
    s[1,0] =  0.04
    s[1,1] = -1e4 * y_n[2] - 6e7 * y_n[1]
    s[1,2] = -1e4 * y_n[1]
    #s[2] =  3e7 * y_n[1]**2
    s[2,1] =  6e7 * y_n[1]
    return s

# Initial condition

def backwards_euler(y0, dt, n_time_steps, n_newton = 5):
    yn = y0.copy()
    ys = [y0]
    time = [0.0]
    I = np.eye(3)
    for t in range(n_time_steps):
        y_guess = yn + 1e-4 * np.ones(3)
        for _ in range(n_newton):
            res = (y_guess - yn) / dt - source_1pt3(y_guess)
            J = I / dt - dsource_1pt3_dy(y_guess)
            dy, info = gmres_custom(J, -res)
            y_guess = y_guess + dy
            if norm2(res) < 1e-10:
                break
        yn = y_guess
        ys.append(yn)
        time.append(time[-1] + dt)

    ys = np.array(ys)
    time = np.array(time)
    return time, ys

def sdirk2(y0, dt, n_time_steps, n_newton=5):
    yn = y0.copy()
    ys = [y0]
    time = [0.0]
    I = np.eye(len(y0))

    gamma = 1.0 - 1.0 / np.sqrt(2)
    one_minus_gamma = 1.0 - gamma

    for t in range(n_time_steps):
        # Stage 1: Solve for k1 in k1 = f(y_n + gamma * dt * k1)
        k1 = source_1pt3(yn)
        for _ in range(n_newton):
            y_stage = yn + gamma * dt * k1
            f_val = source_1pt3(y_stage)
            res = k1 - f_val
            J = I - gamma * dt * dsource_1pt3_dy(y_stage)
            delta, info = gmres_custom(J, -res)
            k1 = k1 + delta
            if np.linalg.norm(res) < 1e-10:
                break

        # Stage 2: Solve for k2 in k2 = f(y_n + (1-gamma)*dt*k1 + gamma*dt*k2)
        k2 = source_1pt3(yn)
        for _ in range(n_newton):
            y_stage = yn + dt * (one_minus_gamma * k1 + gamma * k2)
            f_val = source_1pt3(y_stage)
            res = k2 - f_val
            J = I - gamma * dt * dsource_1pt3_dy(y_stage)
            delta, info = gmres_custom(J, -res)
            k2 = k2 + delta
            if np.linalg.norm(res) < 1e-10:
                break

        yn = yn + dt * (one_minus_gamma * k1 + gamma * k2)
        ys.append(yn)
        time.append(time[-1] + dt)

    return np.array(time), np.array(ys)


def sdirk4(y0, dt, n_time_steps, n_newton=5):
    yn = y0.copy()
    ys = [y0]
    time = [0.0]
    I = np.eye(len(y0))

    # Coefficients for Kaps–Rentrop SDIRK-4 method
    
    #0.20
    #0.15 0.20
    #0.78518518518519 0.98518518518519 0.20
    #0.70671936788942 0.19819311123659 0.09147374364765 0.20
    #0.4074074074074 0.4074074074074 0.10714285714286 0.46851851851852
    
    gamma = .26666666666666666666666666666666670

    a11 = .26666666666666666666666666666666670
    a21 = .50000000000000000000000000000000000
    a22 = .26666666666666666666666666666666670
    a31 = .35415395284327323162274618585298200
    a32 = -.5415395284327323162274618585298197e-1
    a33 = .26666666666666666666666666666666670
    a41 = .8515494131138652076337791881433756e-1
    a42 = -.6484332287891555171683963466229754e-1
    a43 = .7915325296404206392428857585141242e-1
    a44 = .26666666666666666666666666666666670
    a51 = 2.1001157005669327779706120559990740
    a52 = -.76778002844459768133431021850622760
    a53 = 2.3998163610800263980947462052738800
    a54 = -2.9988186998690281613977147094333940
    a55 = .26666666666666666666666666666666670

    b1   = 2.1001157005669327779706120559990740
    b2   = -.76778002844459768133431021850622760
    b3   = 2.3998163610800263980947462052738800
    b4   = -2.9988186998690281613977147094333940
    b5   = .26666666666666666666666666666666670

    #gamma = 0.4358665215  # diagonal elements a_ii for all stages
    #a21 = 0.5529291481
    #a31 = 0.2466725606
    #a32 = 0.4265742951
    #a41 = 0.1881405927
    #a42 = 0.6211338564
    #a43 = -0.0680079704

    #b1  = 0.1881405927
    #b2  = 0.6211338564
    #b3  = -0.0680079704
    #b4  = 0.2587335203

    for t in range(n_time_steps):
        k = [source_1pt3(yn) for _ in range(5)]

        # Stage 1
        for _ in range(n_newton):
            y_stage = yn + dt * gamma * k[0]
            f = source_1pt3(y_stage)
            res = k[0] - f
            J = I - dt * gamma * dsource_1pt3_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[0] += delta
            if norm2(res) < 1e-10:
                break

        # Stage 2
        for _ in range(n_newton):
            y_stage = yn + dt * (a21 * k[0] + gamma * k[1])
            f = source_1pt3(y_stage)
            res = k[1] - f
            J = I - dt * gamma * dsource_1pt3_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[1] += delta
            if norm2(res) < 1e-10:
                break

        # Stage 3
        for _ in range(n_newton):
            y_stage = yn + dt * (a31 * k[0] + a32 * k[1] + gamma * k[2])
            f = source_1pt3(y_stage)
            res = k[2] - f
            J = I - dt * gamma * dsource_1pt3_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[2] += delta
            if norm2(res) < 1e-10:
                break

        # Stage 4
        for _ in range(n_newton):
            y_stage = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2] + gamma * k[3])
            f = source_1pt3(y_stage)
            res = k[3] - f
            J = I - dt * gamma * dsource_1pt3_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[3] += delta
            if norm2(res) < 1e-10:
                break
    
        # Stage 5
        for _ in range(n_newton):
            y_stage = yn + dt * (a51 * k[0] + a52 * k[1] + a53 * k[2] +a54 * k[3] + gamma * k[4])
            f = source_1pt3(y_stage)
            res = k[4] - f
            J = I - dt * gamma * dsource_1pt3_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[4] += delta
            if norm2(res) < 1e-10:
                break
        # Final update
        yn = yn + dt * (b1 * k[0] + b2 * k[1] + b3 * k[2] + b4 * k[3] + b5 * k[4])
        ys.append(yn)
        time.append(time[-1] + dt)

    return np.array(time), np.array(ys)

def rosenbrock2(y0, dt, n_time_steps):
    yn = y0.copy()
    ys = [y0]
    time = [0.0]
    I = np.eye(len(y0))

    gamma = 1.0 + 1.0 / np.sqrt(2)  # ~0.2928932188
    alpha = 1.0/gamma
    beta = -2/gamma
    m1 = 3.0/(2.0*gamma)
    m2 = 1.0/(2*gamma)

    for t in range(n_time_steps):
        J = dsource_1pt3_dy(yn)  # Jacobian at y_n
        G = (1.0/(gamma* dt)) * I - J   # shared for all stages

        # Stage 1
        rhs1 = source_1pt3(yn)
        k1, _ = gmres_custom(G, rhs1)

        # Stage 2
        y_stage = yn + alpha * k1
        rhs2 = source_1pt3(y_stage) + beta/dt *k1 - G @ k1 #source_1pt3(y_stage) - J @ (a21 * dt * k1)
        dk2, _ = gmres_custom(G, rhs2)
        k2 = k1 + dk2

        # Combine stages
        yn = yn + m1 * k1 + m2 * k2
        ys.append(yn)
        time.append(time[-1] + dt)

    return np.array(time), np.array(ys)


y0 = np.array([1.0, 0.0, 0.0])
dt = 1e-3
time_final = 0.3
n_time_steps = int(time_final / dt)

(time_be, y_be), t_be = timeit_solver("Backward Euler", backwards_euler, y0, dt, n_time_steps, n_newton=5)
(time_sdirk2, y_sdirk2), t_sdirk2 = timeit_solver("SDIRK2", sdirk2, y0, dt, n_time_steps, n_newton=5)
(time_sdirk4, y_sdirk4), t_sdirk4 = timeit_solver("SDIRK4", sdirk4, y0, dt, n_time_steps, n_newton=5)
(time_be_h2, y_be_h2), t_be_h2 = timeit_solver("Backward Euler (h/2)", backwards_euler, y0, dt/2.0, 2*n_time_steps, n_newton=5)
(time_sdirk2_h2, y_sdirk2_h2), t_sdirk2_h2 = timeit_solver("SDIRK2 (h/2)", sdirk2, y0, dt/2, 2*n_time_steps, n_newton=5)
(time_ros2, y_ros2), t_ros2 = timeit_solver("Rosenbrock2", rosenbrock2, y0, dt, n_time_steps)
(time_ros2_h2, y_ros2_h2), t_ros2_h2 = timeit_solver("Rosenbrock2", rosenbrock2, y0, dt/2, 2*n_time_steps)
(time_sdirk4_h10, y_sdirk4_h10), t_sdirk4_h10 = timeit_solver("SDIRK4", sdirk4, y0, dt/10, n_time_steps*10, n_newton=5)

plt.plot(time_sdirk4_h10, y_sdirk4_h10[:, 1], '-k', label=f'SDIRK4, dt/10\ntime = {t_sdirk4_h10:.2e}', mfc = "white", lw=4)
plt.plot(time_be, y_be[:, 1], '-or', label=f'Backward Euler\ntime = {t_be:.2e}', mfc = "white")
plt.plot(time_sdirk2, y_sdirk2[:, 1], '-ok', label=f'SDIRK2\ntime = {t_sdirk2:.2e}', mfc = "white")
plt.plot(time_ros2, y_ros2[:, 1], '-og', label=f'Rosenbrock2\ntime = {t_ros2:.2e}', mfc = "white")
plt.plot(time_sdirk4, y_sdirk4[:, 1], '-ob', label=f'SDIRK4\ntime = {t_sdirk4:.2e}', mfc = "white")
plt.plot(time_be_h2, y_be_h2[:, 1], '-^r', label=f'Backward Euler, dt/2\ntime = {t_be_h2:.2e}', mfc = "white")
plt.plot(time_sdirk2_h2, y_sdirk2_h2[:, 1], '-^k', label=f'SDIRK2, dt/2\ntime = {t_sdirk2_h2:.2e}', mfc = "white")
plt.plot(time_ros2_h2, y_ros2_h2[:, 1], '-^g', label=f'Rosenbrock2, dt/2\ntime = {t_ros2_h2:.2e}', mfc = "white")
np.save("time_sdirk4_h10.npy", time_sdirk4_h10)
np.save("y_sdirk4_h10.npy", y_sdirk4_h10)

np.save("time_be.npy", time_be)
np.save("y_be.npy", y_be)

np.save("time_sdirk2.npy", time_sdirk2)
np.save("y_sdirk2.npy", y_sdirk2)

np.save("time_ros2.npy", time_ros2)
np.save("y_ros2.npy", y_ros2)

np.save("time_sdirk4.npy", time_sdirk4)
np.save("y_sdirk4.npy", y_sdirk4)

np.save("time_be_h2.npy", time_be_h2)
np.save("y_be_h2.npy", y_be_h2)

np.save("time_sdirk2_h2.npy", time_sdirk2_h2)
np.save("y_sdirk2_h2.npy", y_sdirk2_h2)

np.save("time_ros2_h2.npy", time_ros2_h2)
np.save("y_ros2_h2.npy", y_ros2_h2)

plt.ylim([3.636e-5, 3.65e-5])
plt.xlim([0.0025, 0.02])
plt.legend(loc=1, ncol=2, fontsize=8)
plt.xlabel("Time")
plt.ylabel("$y_2$")
plt.title("Implicit Solve of dy/dt = S(y)")
plt.grid(True)
plt.show()
