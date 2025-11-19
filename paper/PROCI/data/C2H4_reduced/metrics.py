
#!/usr/bin/env python3
import numpy as np

def departure(A):
    "Count zero rows as a simple non-normality proxy (user's original)."
    return np.sum(~A.any(axis=1))

def inv_largest_tau(J):
    # tau_min = 1 / max(|Re(lambda)|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size == 0 else 1.0 / np.max(nz)

def inv_smallest_tau(J):
    # tau_max = 1 / min(|Re(lambda)|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size == 0 else 1.0 / np.min(nz)

def eigs(J):
    return 1/np.abs(np.real(np.linalg.eigvals(J)))

def svs(A):
    return np.linalg.svd(A, compute_uv=False)

def alphas(A,b):
    #U, s, Vt = np.linalg.svd(A, full_matrices=True)
    #return U.T @ b   
    eigvals, V = np.linalg.eig(A)
    return np.linalg.inv(V) @ b

def stiffness(J):
    eig = np.linalg.eigvals(J)
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size < 2 else np.max(nz) / np.min(nz)

# --- 1-D clustering helpers ---
def nn_distances(x):
    x = np.sort(np.asarray(x, dtype=float))
    if x.size < 2:
        return np.array([np.nan]), x
    d = np.empty_like(x)
    d[0]  = x[1] - x[0]
    d[-1] = x[-1] - x[-2]
    if x.size > 2:
        d[1:-1] = np.minimum(x[1:-1]-x[:-2], x[2:]-x[1:-1])
    return d, x

def clustering_scores(x):
    d, xs = nn_distances(x)
    if xs.size < 2:
        return {"R_NN": np.nan, "CV_NN": np.nan}
    mean_d = np.nanmean(d)
    cv_nn  = np.nanstd(d, ddof=0) / mean_d if mean_d > 0 else np.nan
    norm   = (xs[-1]-xs[0]) / (len(xs)-1) if len(xs) > 1 else np.nan
    R_nn   = mean_d / norm if norm and norm > 0 else np.nan
    return {"R_NN": R_nn, "CV_NN": cv_nn}

# Colorblind-friendly (Okabe–Ito)
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#000000", "#CC79A7", "#F0E442"]
