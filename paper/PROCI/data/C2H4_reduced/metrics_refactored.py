#!/usr/bin/env python3
import numpy as np

def departure(A): return np.sum(~A.any(axis=1))

def inv_largest_tau(J):
    eig = np.linalg.eigvals(J); mag = np.abs(np.real(eig))
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]; return np.nan if nz.size == 0 else 1.0 / np.max(nz)

def inv_smallest_tau(J):
    eig = np.linalg.eigvals(J); 
    mag = np.abs(np.real(eig))
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]; return np.nan if nz.size == 0 else 1.0 / np.min(nz)

def eigs(J):
    return np.linalg.eigvals(J)

def svs(A): return np.linalg.svd(A, compute_uv=False)

import numpy as np
def alphas(A,b):
    #dead
    return b
def contribution_to_b(A, b, k=5, j=5):
    eigvals, V = np.linalg.eig(A)
    #V = V / np.linalg.norm(V, axis=0, keepdims=True)

    # Coefficients of b in eigenbasis
    c = np.linalg.solve(V, b)
    w = np.abs(c)
    #w /= np.sum(w)

    # Sort by eigenvalue magnitude
    idx = np.argsort(np.abs(eigvals))
    w = w[idx]

    n = len(w)
    k = min(k, n)
    j = min(j, n)


    return np.abs(w)#np.array([float(bottom_frac), float(top_frac)])

from scipy.linalg import eig  # preferred: returns left & right
def modal_contribution(A, b, sort='abs', k=5, j=5):

    eigvals, VL, VR = eig(A, left=True, right=True)
    # Biorthonormalize so u_i^H v_i = 1
    for i in range(VR.shape[1]):
        s = VL[:, i].conj().T @ VR[:, i]
        VR[:, i] /= s  # now VL[:,i]^H VR[:,i] = 1

    # Modal projections
    ub = VL.conj().T @ b            # u_i^H b
    Pib = VR * ub[np.newaxis, :]    # each column i is v_i * (u_i^H b)
    w = np.sum(np.abs(Pib)**2, axis=0)  # ||P_i b||^2
    w /= (np.linalg.norm(b)**2 + 0.0)
    w /=np.sum(w)

    # sort
    if sort == 'abs':
        idx = np.argsort(np.abs(eigvals))
    elif sort == 'real':
        idx = np.argsort(np.real(eigvals))
    else:
        idx = np.arange(len(eigvals))
    return  w[idx] #eigvals[idx],

import numpy as np
import scipy.linalg as la

import numpy as np
import scipy.linalg as la

def schur_alignment(A, b, k=5, j=5, sort_by="real"):
    # Schur form: A = Q T Q^H (Q unitary)
    scale = np.linalg.norm(A, ord=np.inf)
    A = A/scale
    b = b/scale
    T, Q = la.schur(A, output='complex')
    lam = np.diag(T)

    if sort_by == "abs":
        order = np.argsort(np.abs(lam))
    elif sort_by == "real":
        order = np.argsort(np.real(lam))
    elif sort_by == "imag":
        order = np.argsort(np.imag(lam))
    else:
        raise ValueError("sort_by must be 'abs', 'real', or 'imag'.")

    Qs   = Q[:, order]
    # Coordinates of b in the Schur basis
    bhat = Qs.conj().T @ b
    w    = np.abs(bhat)**2             # energy per Schur vector

    denom = np.linalg.norm(b)**2
    if denom == 0.0:
        return np.array([0.0, 0.0])

    w /= denom                         # sum(w) ≈ 1

    n = len(w)
    k = max(0, min(int(k), n))
    j = max(0, min(int(j), n))

    bottom_frac = float(np.sum(w[:k]))       # smallest by chosen ordering
    top_frac    = float(np.sum(w[n-j:]))     # largest by chosen ordering
    return np.array([bottom_frac, top_frac])

def liao_metric(A, b, k=5, j=5):
    eigvals, V = np.linalg.eig(A)
    a = np.linalg.inv(V) @ b
    C = np.zeros_like(A, dtype=complex)
    for i in range(len(a)):
        C[:,i] = a[i] * V[:,i]
    svds = np.linalg.svd(C, compute_uv=False)
    return np.array([np.max(svds), 0])

def svd_alignment(A, b, k=5, j=5):
    U, _, _ = np.linalg.svd(A)
    a = U.T @ b
    return a

def stiffness(J):
    eig = np.linalg.eigvals(J); mag = np.abs(np.real(eig))
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]; return np.nan if nz.size < 2 else np.max(nz) / np.min(nz)

PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#000000", "#CC79A7", "#F0E442"]
