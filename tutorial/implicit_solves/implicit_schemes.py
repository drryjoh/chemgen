import numpy as np
from linear_solve import gmres_custom
from source_functions import source, dsource_dy
from my_math import norm2, least_squares_qr, back_substitution

def backwards_euler(y0, dt, n_time_steps, n_newton = 5):
    yn = y0.copy()
    ys = [y0]
    time = [0.0]
    I = np.eye(3)
    for t in range(n_time_steps):
        y_guess = yn + 1e-4 * np.ones(3)
        for _ in range(n_newton):
            res = (y_guess - yn) / dt - source(y_guess)
            J = I / dt - dsource_dy(y_guess)
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
        k1 = source(yn)
        for _ in range(n_newton):
            y_stage = yn + gamma * dt * k1
            f_val = source(y_stage)
            res = k1 - f_val
            J = I - gamma * dt * dsource_dy(y_stage)
            delta, info = gmres_custom(J, -res)
            k1 = k1 + delta
            if np.linalg.norm(res) < 1e-10:
                break

        # Stage 2: Solve for k2 in k2 = f(y_n + (1-gamma)*dt*k1 + gamma*dt*k2)
        k2 = source(yn)
        for _ in range(n_newton):
            y_stage = yn + dt * (one_minus_gamma * k1 + gamma * k2)
            f_val = source(y_stage)
            res = k2 - f_val
            J = I - gamma * dt * dsource_dy(y_stage)
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

    for t in range(n_time_steps):
        k = [source(yn) for _ in range(5)]

        # Stage 1
        for _ in range(n_newton):
            y_stage = yn + dt * gamma * k[0]
            f = source(y_stage)
            res = k[0] - f
            J = I - dt * gamma * dsource_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[0] += delta
            if norm2(res) < 1e-10:
                break

        # Stage 2
        for _ in range(n_newton):
            y_stage = yn + dt * (a21 * k[0] + gamma * k[1])
            f = source(y_stage)
            res = k[1] - f
            J = I - dt * gamma * dsource_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[1] += delta
            if norm2(res) < 1e-10:
                break

        # Stage 3
        for _ in range(n_newton):
            y_stage = yn + dt * (a31 * k[0] + a32 * k[1] + gamma * k[2])
            f = source(y_stage)
            res = k[2] - f
            J = I - dt * gamma * dsource_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[2] += delta
            if norm2(res) < 1e-10:
                break

        # Stage 4
        for _ in range(n_newton):
            y_stage = yn + dt * (a41 * k[0] + a42 * k[1] + a43 * k[2] + gamma * k[3])
            f = source(y_stage)
            res = k[3] - f
            J = I - dt * gamma * dsource_dy(y_stage)
            delta, _ = gmres_custom(J, -res)
            k[3] += delta
            if norm2(res) < 1e-10:
                break
    
        # Stage 5
        for _ in range(n_newton):
            y_stage = yn + dt * (a51 * k[0] + a52 * k[1] + a53 * k[2] +a54 * k[3] + gamma * k[4])
            f = source(y_stage)
            res = k[4] - f
            J = I - dt * gamma * dsource_dy(y_stage)
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
        J = dsource_dy(yn)  # Jacobian at y_n
        G = (1.0/(gamma* dt)) * I - J   # shared for all stages

        # Stage 1
        rhs1 = source(yn)
        k1, _ = gmres_custom(G, rhs1)

        # Stage 2
        y_stage = yn + alpha * k1
        rhs2 = source(y_stage) + beta/dt *k1 - G @ k1 #source(y_stage) - J @ (a21 * dt * k1)
        dk2, _ = gmres_custom(G, rhs2)
        k2 = k1 + dk2

        # Combine stages
        yn = yn + m1 * k1 + m2 * k2
        ys.append(yn)
        time.append(time[-1] + dt)

    return np.array(time), np.array(ys)

def seulex_adaptive(y0, dt_init, t_final, rtol=1e-6, atol=1e-10, n_stages=3, n_newton=5):
    yn = y0.copy()
    ys = [y0]
    time = [0.0]
    I = np.eye(len(y0))

    t = 0.0
    dt = dt_init

    while t < t_final:
        if t + dt > t_final:
            dt = t_final - t

        # Extrapolation table
        U = [np.copy(yn) for _ in range(n_stages)]
        m_vals = [k + 1 for k in range(n_stages)]

        # Stage loop
        for k, m in enumerate(m_vals):
            h = dt / m
            y_k = yn.copy()

            for _ in range(m):
                y_guess = y_k + 1e-4 * np.ones_like(y_k)
                for _ in range(n_newton):
                    res = (y_guess - y_k) / h - source(y_guess)
                    J = I / h - dsource_dy(y_guess)
                    dy, _ = gmres_custom(J, -res)
                    y_guess += dy
                    if np.linalg.norm(res) < 1e-10:
                        break
                y_k = y_guess

            U[k] = y_k

        # Extrapolation (Neville–Aitken)
        for j in range(1, n_stages):
            for k in range(n_stages - 1, j - 1, -1):
                factor = m_vals[k] / m_vals[k - j] - 1.0
                U[k] = U[k] + (U[k] - U[k - 1]) / factor

        y_extrapolated = U[n_stages - 1]
        y_lower = U[n_stages - 2]

        # Error estimation
        scale = atol + rtol * np.maximum(np.abs(y_extrapolated), np.abs(yn))
        error = np.linalg.norm((y_extrapolated - y_lower) / scale) / np.sqrt(len(y0))

        # Accept step
        if error <= 1.0:
            t += dt
            yn = y_extrapolated
            ys.append(yn)
            time.append(t)

        # Compute next dt (Hairer's suggestion)
        safety = 0.9
        if error == 0.0:
            dt = dt * 2.0
        else:
            dt = dt * min(5.0, max(0.1, safety * error ** (-1.0 / (n_stages + 1))))

    return np.array(time), np.array(ys)

