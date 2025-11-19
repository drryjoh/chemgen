#!/usr/bin/env python3
import math
import numpy as np
from scipy.sparse.linalg import gmres, lsqr, LinearOperator
import eigen_gmres

def error_norm(dx, x_ref, abs_tol, rel_tol):
    denom = abs_tol + rel_tol * np.abs(x_ref)
    return np.linalg.norm(dx / denom)

def full_jacobian(J, dt):
    I = np.eye(J.shape[0])
    A = I / dt - J
    return A

def full_jacobian_scaled(J, dt):
    I = np.eye(J.shape[0])
    A = I / dt - J
    A_scaled = A 
    scale = np.linalg.norm(A, ord=np.inf)
    A_scaled = A / scale
    return A_scaled

def _diag_inv(A):
    diag = np.diag(A)
    with np.errstate(divide="ignore", invalid="ignore"):
        invd = np.where(np.abs(diag) > 0, 1.0 / diag, 1.0)
    return invd

def gmres_with_iterations(A, b, tol=1e-4, atol=1e-7, maxiter=100):
    M_inv = _diag_inv(A)
    def Mv(v): return M_inv * v
    M = LinearOperator(A.shape, matvec=Mv)
    iters = [0]
    def count(_): iters[0] += 1
    x, info = gmres(A, b, M=M, rtol=tol, atol=atol, maxiter=maxiter, callback=count)
    return x, info, iters[0]

def lsqr_with_iterations(A, b, tol=1e-4, maxiter=100, use_preconditioner=False):
    A = np.asarray(A); b = np.asarray(b).ravel()
    if use_preconditioner:
        M_inv = _diag_inv(A)
        def Mv(v): return M_inv * v
        def AMv(v): return A @ Mv(v)
        def AMTv(w): return Mv(A.T @ w)
        A_pre = LinearOperator(A.shape, matvec=AMv, rmatvec=AMTv)
        result = lsqr(A_pre, b, atol=tol, btol=tol, iter_lim=maxiter, show=False)
        y, istop, itn = result[0], result[1], result[2]
        x = Mv(y)
    else:
        result = lsqr(A, b, atol=tol, btol=tol, iter_lim=maxiter, show=False)
        x, istop, itn = result[0], result[1], result[2]
    return x, istop, itn
def sdirk2(C, T, dt, be, gmres_tolerance=1e-6, max_iter=10, gmres_method="numpy"):
    abs_newton_tol, rel_newton_tol = 1e-7, 1e-6
    mod = be["mod"]
    y = np.zeros(len(C) + 1); y[1:] = C
    y[0] = mod.internal_energy_volume_specific(C, T) if be.get("require_internal_energy", "no") == "yes" else T
    y_init = y.copy()
    gamma = 1.0 - 1.0 / math.sqrt(2.0); one_minus_gamma = 1.0 - gamma
    T0 = be["temp_from_state"](y); k1 = np.array(mod.source(y_init[1:], T0))
    I = np.eye(len(k1))
    running_newton = linear_iterations = 0

    for _ in range(max_iter):
        y_stage = y_init + gamma * dt * np.concatenate(([0.0], k1))
        T_ = be["temp_from_state"](y_stage)
        f_val = np.array(mod.source(y_stage[1:], T_))
        res = f_val - k1
        J = I - gamma * dt * np.array(mod.source_jacobian(y_stage[1:], T_))
        solver = gmres_with_iterations if gmres_method == "numpy" else lsqr_with_iterations
        dk, _, nlin = solver(J, res, tol=gmres_tolerance, maxiter=100)
        k1 += dk; linear_iterations += nlin; running_newton += 1
        if error_norm(dk, k1, abs_newton_tol, rel_newton_tol) < 1.0: break

    k2 = k1.copy()
    for _ in range(max_iter):
        y_stage = y_init + dt * (one_minus_gamma * np.concatenate(([0.0], k1)) + gamma * np.concatenate(([0.0], k2)))
        T_ = be["temp_from_state"](y_stage)
        f_val = np.array(mod.source(y_stage[1:], T_))
        res = f_val - k2
        J = I - gamma * dt * np.array(mod.source_jacobian(y_stage[1:], T_))
        dk, _, nlin = gmres_with_iterations(J, res, tol=gmres_tolerance, maxiter=100)
        k2 += dk; linear_iterations += nlin; running_newton += 1
        if error_norm(dk, k2, abs_newton_tol, rel_newton_tol) < 1.0: break

    y_next = y_init + dt * (one_minus_gamma * np.concatenate(([0.0], k1)) + gamma * np.concatenate(([0.0], k2)))
    return y_next, linear_iterations, running_newton

def scale_system_A(A, be):
    #n = A.shape[0]
    #scale = np.ones(n, dtype=float)
    #row_scale = np.ones(n, dtype=float)

    #if be["key"] == "cgt":
    #    # column scaling (normalization)
    #    scale[0]  = 1000.0     # divide top column/row by this
    #    scale[1:] = 0.002

    #    # row scaling (pre-multiply rows)
    #    row_scale[0]  = 1000.0  # multiply first row by 1000
    #    row_scale[1:] = 0.002   # multiply remaining rows by 0.002

    #else:  # cgc
    #    scale[0]  = 1.0
    #    scale[1:] = 0.002

    #    row_scale[0]  = 1.0
    #    row_scale[1:] = 0.002

    ## Apply row scaling (multiply) and column normalization (divide)
    #A_scaled = (A * row_scale[:, None]) / scale[None, :]
    return A


def scale_system_b(b, be):
    #if be["key"] == "cgt":
    #    # --- Case 1: top element divided by 1000, rest by 0.002
    #    scale = np.ones_like(b)
    #    scale[0] = 1000.0
    #    scale[1:] = 0.002
    #    b_scaled = b / scale
    #else:
    #    # --- Case 2: top element unchanged, rest divided by 0.002
    #    scale = np.ones_like(b)
    #    scale[1:] = 0.002
    #    b_scaled = b / scale
    return b

def scale_back_b(b, be):
    #if be["key"] == "cgt":
    #    # --- Case 1: top element divided by 1000, rest by 0.002
    #    scale = np.ones_like(b)
    #    scale[0] = 1000.0
    #    scale[1:] = 0.002
    #    b_scaled = b * scale
    #else:
    #    # --- Case 2: top element unchanged, rest divided by 0.002
    #    scale = np.ones_like(b)
    #    scale[1:] = 0.002
    #    b_scaled = b * scale
    return b

def backwards_euler(C, T, dt, be, gmres_tolerance=1e-6, max_iter=10, gmres_method="numpy", linear_solver_verbose = False, step = 0, formulation = ""):
    abs_newton_tol, rel_newton_tol = 1e-7, 1e-6
    mod = be["mod"]
    y = np.zeros(len(C) + 1); y[1:] = C
    y[0] = mod.internal_energy_volume_specific(C, T) if be.get("require_internal_energy", "no") == "yes" else T
    y_init = y.copy(); y_guess = y.copy(); I = np.eye(y.size)
    running_newton = linear_iterations = 0
    if linear_solver_verbose:
        linear_history = []


    for _ in range(max_iter):
        T_ = be["temp_from_state"](y_guess)
        f = np.array(mod.source(y_guess[1:], T_))
        J = np.array(mod.source_jacobian(y_guess[1:], T_))
        A = I / dt -scale_system_A(J,be); res = scale_system_b(f, be) - (1.0 / dt) * scale_system_b((y_guess - y_init), be)
        solver = gmres_with_iterations if gmres_method == "numpy" else lsqr_with_iterations

        A_scaled = A 
        res_scaled = res
        dy, _, it_lin = solver(A_scaled, res_scaled, tol=gmres_tolerance, atol = 1e-12, maxiter=100)
        if linear_solver_verbose:
            np.save(f"data/A_b/{formulation}_{step}_{running_newton}_A.npy", A)
            np.save(f"data/A_b/{formulation}_{step}_{running_newton}_b.npy", res)
            linear_history.append(it_lin)
        y_guess += dy; linear_iterations += it_lin; running_newton += 1
        if error_norm(dy, y_guess, abs_newton_tol, rel_newton_tol) < 1.0: break
    if linear_solver_verbose:
        np.save(f"data/A_b/{formulation}_{step}_linear_history.npy", np.array(linear_history))
    return y_guess, linear_iterations, running_newton
