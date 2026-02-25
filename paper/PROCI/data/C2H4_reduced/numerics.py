
#!/usr/bin/env python3
import math
import numpy as np
from scipy.sparse.linalg import gmres, lsqr, LinearOperator
import eigen_gmres

def error_norm(x, x_ref, abs_tolerance, rel_tolerance):
    denom = abs_tolerance + rel_tolerance * np.abs(x_ref)
    return np.linalg.norm(x / denom)

def full_jacobian(J, dt):
    # Backward Euler form kept as identity on J for now
    I = np.eye(np.shape(J)[0])
    return I/dt - J

def jacobi_precondition(A, eps=1e-14):
    """Return D^{-1}A and D^{-1} itself."""
    d = np.diag(A)
    nz = np.abs(d) > eps
    invd = np.zeros_like(d)
    invd[nz] = 1.0 / d[nz]
    M_inv = np.diag(invd)
    return M_inv @ A

def gmres_with_iterations(A, b, tol=1e-4, maxiter=100):
    """Right-preconditioned GMRES with a diagonal preconditioner and iteration counter."""
    diag_A = np.diag(A)
    # avoid division-by-zero; fall back to ones where needed
    with np.errstate(divide="ignore", invalid="ignore"):
        M_inv_vec = np.where(np.abs(diag_A) > 0, 1.0 / diag_A, 1.0)

    def Mv(v):
        return M_inv_vec * v

    M = LinearOperator(A.shape, matvec=Mv)
    iters = [0]
    def counter(residual):
        iters[0] += 1

    x, info = gmres(A, b, M=M, rtol=tol, atol=1e-7, maxiter=maxiter, callback=counter)
    return x, info, iters[0]

def lsqr_with_iterations(A, b, tol=1e-4, maxiter=100, use_preconditioner=False):
    """
    LSQR solver with optional right diagonal preconditioning and iteration counter.

    Parameters
    ----------
    A : array_like or LinearOperator
        System matrix.
    b : array_like
        Right-hand side vector.
    tol : float, optional
        Convergence tolerance for LSQR (default 1e-4).
    maxiter : int, optional
        Maximum number of iterations (default 100).
    use_preconditioner : bool, optional
        If True, applies right Jacobi preconditioning; if False, solves unpreconditioned system.
    """
    A = np.asarray(A)
    b = np.asarray(b).ravel()

    if use_preconditioner:
        # Jacobi (diagonal) right preconditioner
        diag_A = np.diag(A)
        with np.errstate(divide="ignore", invalid="ignore"):
            M_inv_vec = np.where(np.abs(diag_A) > 0, 1.0 / diag_A, 1.0)

        def Mv(v):
            return M_inv_vec * v

        M = LinearOperator(A.shape, matvec=Mv, rmatvec=Mv)

        # Define preconditioned operator: A_pre = A * M
        def AMv(v):
            return A @ Mv(v)

        def AMTv(w):
            return Mv(A.T @ w)

        A_pre = LinearOperator(A.shape, matvec=AMv, rmatvec=AMTv)

        # Solve (A * M) y = b
        result = lsqr(A_pre, b, atol=tol, btol=tol, iter_lim=maxiter, show=False)
        y, istop, itn = result[0], result[1], result[2]

        # Recover x = M * y
        x = Mv(y)
    else:
        # Unpreconditioned solve
        result = lsqr(A, b, atol=tol, btol=tol, iter_lim=maxiter, show=False)
        x, istop, itn = result[0], result[1], result[2]

    return x, istop, itn


def sdirk2(C, T, dt, be, gmres_tolerance=1e-6, max_iter=10, gmres_method="numpy"):
    """Two-stage SDIRK(2) step. Returns (y_next, total_linear_iters, newton_iters)."""
    abs_newton_tol = 1e-7
    rel_newton_tol = 1e-6

    mod = be["mod"]
    y = np.zeros(len(C) + 1)
    y[1:] = C
    if be.get("require_internal_energy", "no") == "yes":
        y[0] = float(mod.internal_energy_volume_specific(C, T))
    else:
        y[0] = float(T)

    y_init = y.copy()
    gamma = 1.0 - 1.0 / math.sqrt(2.0)
    one_minus_gamma = 1.0 - gamma

    T0 = be["temp_from_state"](y)
    k1 = np.array(mod.source(y_init[1:], T0), dtype=float)
    I = np.eye(len(k1), dtype=float)

    running_newton = 0
    linear_iterations = 0

    # Stage 1
    for _ in range(max_iter):
        y_stage = y_init + gamma * dt * np.concatenate(([0.0], k1))
        T_ = be["temp_from_state"](y_stage)
        f_val = np.array(mod.source(y_stage[1:], T_), dtype=float)
        res = f_val - k1
        J = I - gamma * dt * np.array(mod.source_jacobian(y_stage[1:], T_), dtype=float)

        if gmres_method == "numpy":
            dk, info, nlin = gmres_with_iterations(J, res, tol=gmres_tolerance, maxiter=100)
        elif gmres_method == "lsqr":
            dk, info, nlin = lsqr_with_iterations(J, res, tol=gmres_tolerance, maxiter=100)
        else:
            r = eigen_gmres.gmres_dense(J, res, rtol=gmres_tolerance, maxiter=100, restart=10)
            dk, info, nlin = r.x, r.info, r.iters

        linear_iterations += nlin
        if info > 0:
            print(f"SDIRK stage-1 GMRES did not converge after {info} iterations")
        k1 = k1 + dk
        running_newton += 1
        if error_norm(dk, k1, abs_newton_tol, rel_newton_tol) < 1.0:
            break

    # Stage 2
    k2 = k1.copy()
    for _ in range(max_iter):
        y_stage = y_init + one_minus_gamma * dt * np.concatenate(([0.0], k1)) + gamma * dt * np.concatenate(([0.0], k2))
        T_ = be["temp_from_state"](y_stage)
        f_val = np.array(mod.source(y_stage[1:], T_), dtype=float)
        res = f_val - k2
        J = I - gamma * dt * np.array(mod.source_jacobian(y_stage[1:], T_), dtype=float)

        if gmres_method == "numpy":
            dk, info, nlin = gmres_with_iterations(J, res, tol=gmres_tolerance, maxiter=100)
        else:
            r = eigen_gmres.gmres_dense(J, res, rtol=gmres_tolerance, maxiter=100, restart=10)
            dk, info, nlin = r.x, r.info, r.iters

        linear_iterations += nlin
        if info > 0:
            print(f"SDIRK stage-2 GMRES did not converge after {info} iterations")
        k2 = k2 + dk
        running_newton += 1
        if error_norm(dk, k2, abs_newton_tol, rel_newton_tol) < 1.0:
            break

    y_next = y_init + dt * (one_minus_gamma * np.concatenate(([0.0], k1)) + gamma * np.concatenate(([0.0], k2)))
    return y_next, linear_iterations, running_newton

def backwards_euler(C, T, dt, be, gmres_tolerance=1e-6, max_iter=10, gmres_method="numpy"):
    """
    Implicit backward-Euler step:
      Solve [(I/dt) - J(y_guess)] * dy = f(y_guess) - (1/dt) * (y_guess - y_init)
      y_{n+1} = y_guess + dy, iterate Newton until error_norm < 1.
    Returns:
      y_next, total_linear_iters, newton_iters
    """
    abs_newton_tol = 1e-7
    rel_newton_tol = 1e-6

    mod = be["mod"]
    y = np.zeros(len(C) + 1, dtype=float)
    y[1:] = C
    if be.get("require_internal_energy", "no") == "yes":
        y[0] = float(mod.internal_energy_volume_specific(C, T))
    else:
        y[0] = float(T)

    y_init = y.copy()
    y_guess = y.copy()

    nvars = y_init.size
    I = np.eye(nvars, dtype=float)

    running_newton = 0
    linear_iterations = 0

    for _ in range(max_iter):
        T_ = be["temp_from_state"](y_guess)
        f = np.array(mod.source(y_guess[1:], T_), dtype=float)
        J = np.array(mod.source_jacobian(y_guess[1:], T_), dtype=float)
        A = I / dt - J
        res = f - (1.0 / dt) * (y_guess - y_init)

        # Solve
        if gmres_method == "numpy":
            dy_sys, info, iters_lin = gmres_with_iterations(A, res, tol=gmres_tolerance, maxiter=100)
        elif gmres_method == "lsqr":
            dy_sys, info, iters_lin = lsqr_with_iterations(A, res, tol=gmres_tolerance/100, maxiter=200)
        else:
            r = eigen_gmres.gmres_dense(A, res, rtol=gmres_tolerance, maxiter=100, restart=50)
            dy_sys, info, iters_lin = r.x, r.info, r.iters


        dy = dy_sys

        linear_iterations += iters_lin
        if info > 0:
            print(f"GMRES did not converge after {linear_iterations} iterations: info={info}")

        y_guess = y_guess + dy
        running_newton += 1

        if error_norm(dy, y_guess, abs_newton_tol, rel_newton_tol) < 1.0:
            break

    return y_guess, linear_iterations, running_newton
