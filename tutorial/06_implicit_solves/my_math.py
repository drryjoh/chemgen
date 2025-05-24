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

