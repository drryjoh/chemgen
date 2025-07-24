import numpy as np
from my_math import back_substitution, least_squares_qr, norm2

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

