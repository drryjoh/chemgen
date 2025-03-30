#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import gmres

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
            dy, info = gmres(J, -res)
            y_guess = y_guess + dy
            if np.linalg.norm(res) < 1e-10:
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

    gamma = 1.0 / np.sqrt(2)
    one_minus_gamma = 1.0 - gamma

    for t in range(n_time_steps):
        # Stage 1: Solve for k1 in k1 = f(y_n + gamma * dt * k1)
        k1 = source_1pt3(yn)
        for _ in range(n_newton):
            y_stage = yn + gamma * dt * k1
            f_val = source_1pt3(y_stage)
            res = k1 - f_val
            J = I - gamma * dt * dsource_1pt3_dy(y_stage)
            delta, info = gmres(J, -res)
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
            delta, info = gmres(J, -res)
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

    gamma = 0.4358665215

    a21 = 0.5529291481
    a31 = 0.2466725606
    a32 = 0.4265742951
    a41 = 0.1881405927
    a42 = 0.6211338564
    a43 = -0.0680079704
    b1  = 0.1881405927
    b2  = 0.6211338564
    b3  = -0.0680079704
    b4  = 0.2587335203

    for t in range(n_time_steps):
        k = [np.zeros_like(yn) for _ in range(4)]

        # Stage 1
        k[0] = source_1pt3(yn)
        for _ in range(n_newton):
            y_stage = yn + gamma * dt * k[0]
            res = k[0] - source_1pt3(y_stage)
            J = I - gamma * dt * dsource_1pt3_dy(y_stage)
            delta, _ = gmres(J, -res)
            k[0] += delta
            if np.linalg.norm(res) < 1e-10:
                break

        # Stage 2
        y2 = yn + dt * (a21 * k[0])
        k[1] = source_1pt3(y2)
        for _ in range(n_newton):
            y_stage = yn + dt * (a21 * k[0] + gamma * k[1])
            res = k[1] - source_1pt3(y_stage)
            J = I - gamma * dt * dsource_1pt3_dy(y_stage)
            delta, _ = gmres(J, -res)
            k[1] += delta
            if np.linalg.norm(res) < 1e-10:
                break

        # Stage 3
        y3 = yn + dt * (a31 * k[0] + a32 * k[1])
        k[2] = source_1pt3(y3)
        for _ in range(n_newton):
            y_stage = yn + dt * (a31 * k[0] + a32 * k[1] + gamma * k[2])
            res = k[2] - source_1pt3(y_stage)
            J = I - gamma * dt * dsource_1pt3_dy(y_stage)
            delta, _ = gmres(J, -res)
            k[2] += delta
            if np.linalg.norm(res) < 1e-10:
                break

        # Stage 4
        y4 = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2])
        k[3] = source_1pt3(y4)
        for _ in range(n_newton):
            y_stage = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2] + gamma * k[3])
            res = k[3] - source_1pt3(y_stage)
            J = I - gamma * dt * dsource_1pt3_dy(y_stage)
            delta, _ = gmres(J, -res)
            k[3] += delta
            if np.linalg.norm(res) < 1e-10:
                break

        # Final update
        yn = yn + dt * (b1 * k[0] + b2 * k[1] + b3 * k[2] + b4 * k[3])
        ys.append(yn)
        time.append(time[-1] + dt)

    return np.array(time), np.array(ys)

y0 = np.array([1.0, 0.0, 0.0])
dt = 1e-2
time_final = 0.3
n_time_steps = int(time_final / dt)

#time_be, y_be = backwards_euler(y0, dt, n_time_steps, n_newton = 5)
#time_be_h2, y_be_h2 = backwards_euler(y0, dt/2.0, 2 * n_time_steps, n_newton = 5)
#time_be_h4, y_be_h4 = backwards_euler(y0, dt/4.0, 4 * n_time_steps, n_newton = 5)
time_be_h8, y_be_h8 = backwards_euler(y0, dt/8.0, 8 * n_time_steps, n_newton = 5)
time_sdirk2, y_sdirk2 = sdirk2(y0, dt/8.0, 8*n_time_steps, n_newton = 5)
time_sdirk4, y_sdirk4 = sdirk4(y0, dt/8.0, 8*n_time_steps, n_newton = 5)
time_be_h16, y_be_h16 = backwards_euler(y0, dt/16.0, 16 * n_time_steps, n_newton = 5)
time_be_h32, y_be_h32 = backwards_euler(y0, dt/100.0, 100 * n_time_steps, n_newton = 5)
time_sdirk2_16, y_sdirk2_16 = sdirk2(y0, dt/16.0, 16*n_time_steps, n_newton = 5)
# Plot y[1] (middle component)
#plt.plot(time_be, y_be[:, 1])
#plt.plot(time_be_h2, y_be_h2[:, 1])
#plt.plot(time_be_h4, y_be_h4[:, 1])
plt.plot(time_be_h8, y_be_h8[:, 1],'-r')
plt.plot(time_sdirk2, y_sdirk2[:, 1],'-k')
plt.plot(time_sdirk4, y_sdirk4[:, 1],'-b')
plt.plot(time_be_h16, y_be_h16[:, 1],'--r')
plt.plot(time_be_h32, y_be_h32[:, 1],'-.r')
plt.plot(time_sdirk2_16, y_sdirk2_16[:, 1],'--k')
plt.ylim([3.625e-5, 3.66e-5])
plt.xlim([0, 0.05])
plt.xlabel("Time")
plt.ylabel("y[1]")
plt.title("Implicit GMRES Solve of dy/dt = source(y)")
plt.grid(True)
plt.show()
