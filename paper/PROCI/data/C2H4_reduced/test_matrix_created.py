#!python3
# Construct A and b from prescribed singular values and left-singular weights (alpha),
# then solve Ax=b using SciPy's GMRES (with a pure NumPy fallback if SciPy is unavailable).
import numpy as np
from scipy.sparse.linalg import gmres as sp_gmres

def random_orthogonal(n, rng):
    M = rng.normal(size=(n, n))
    Q, _ = np.linalg.qr(M)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q

# Problem specification
s = np.array([1e6, 1e3, 1, 1e6], dtype=float)                  # singular values
alpha = np.array([1e6, 1e-6, 1e-6, 1e6], dtype=float)             # left-singular weights
Sigma = np.diag(s)

rng = np.random.default_rng(0)

# Case 1: normal/SPD-like (U=V)
Q = random_orthogonal(4, rng)
U = V = Q
A_spd = U @ Sigma @ V.T
b_spd = U @ alpha
x_true_spd = V @ (alpha / s)

# Case 2: non-normal (U != V)
U2 = random_orthogonal(4, rng)
V2 = random_orthogonal(4, rng)
A_nn = U2 @ Sigma @ V2.T
b_nn = U2 @ alpha
x_true_nn = V2 @ (alpha / s)

def run_case(A, b, x_true, name):
    print(f"--- {name} ---")
    res_hist = []
    def cb(residual):
        res_hist.append(float(residual))
    x, info = sp_gmres(A, b, rtol=1e-12, atol=0.0, callback=cb, maxiter=2000)
    x = np.asarray(x)
    if len(res_hist) == 0:
        res_hist = [np.linalg.norm(b - A @ x)]
    print("SciPy GMRES used.")
    rnorm = np.linalg.norm(b - A @ x)
    err = np.linalg.norm(x - x_true)
    print(f"iterations: {len(res_hist)}")
    print(f"final residual norm: {rnorm:.3e}")
    print(f"solution error norm: {err:.3e}")
    if len(res_hist) > 0:
        print(f"first residual: {res_hist[0]:.3e}, last residual: {res_hist[-1]:.3e}")
    print()

# Run both cases
run_case(A_spd, b_spd, x_true_spd, "Case 1: normal (U=V)")
run_case(A_nn, b_nn, x_true_nn, "Case 2: non-normal (U!=V)")

