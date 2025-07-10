import numpy as np
from scipy.linalg import lu
from scipy import sparse
from scipy.sparse.linalg import splu, spilu
# from splu_helpers import *
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
folder = SCRIPT_DIR.split("chemgen")[0] + "chemgen/bin/modules/"
sys.path.append(os.path.dirname(folder))

from splu_helpers import *

eps = np.finfo(float).eps

def assert_close(A, B, tol=1e-10):
    # absolute tolerance
    try:
        A = A.toarray()
    except AttributeError:
        pass
    try:
        B = B.toarray()
    except AttributeError:
        pass

    if np.linalg.norm(A - B) > tol:
        raise ValueError("np.linalg.norm(A - B) = {}".format(np.linalg.norm(A - B)))

def lu_and_sparsity(A):
    # assert(A.dtype == int or A.dtype == np.int32)

    n = A.shape[0]
    assert(n == A.shape[1])

    sp = np.zeros_like(A, dtype=int) # sparsity
    LU = A.astype(float) # Store L and U within LU

    # First column is not traversed in below for loop (A is copied to LU)
    sp[:,0] = A[:,0] != 0

    def is_nonzero(sp_ij, LU_ij, LU_ik, LU_kj):
        return (sp_ij != 0 or LU_ij != 0 or (LU_ik != 0 and LU_kj != 0))

    for j in range(n):
        # U
        for i in range(j):
            for k in range(i):
                sp[i][j] = is_nonzero(sp[i][j], LU[i][j], LU[i][k], LU[k][j])
                LU[i][j] = LU[i][j] - (LU[i][k] * LU[k][j])

        pivot = -1

        # diag, L
        for i in range(j,n):
            for k in range(j):
                sp[i][j] = is_nonzero(sp[i][j], LU[i][j], LU[i][k], LU[k][j])
                LU[i][j] = LU[i][j] - (LU[i][k] * LU[k][j])

            # print("LU[{}][{}] = {}".format(i, j, LU[i][j]))
            if pivot == -1 and LU[i][j] != 0:
                pivot = i

        if pivot == -1:
            raise ValueError("pivot = {}".format(pivot))
        elif pivot != j:
            raise ValueError("pivot = {}".format(pivot))

        # Only allow diagonal pivots
        if LU[j][j] == 0:
            raise ValueError("j = {}, LU[j][j] = {}".format(j, LU[j][j]))
        if sp[j][j] == 0:
            raise ValueError("j = {}, sp[j][j] = {}".format(j, sp[j][j]))

        # L
        scale = 1./LU[j][j]
        for i in range(j+1,n):
            LU[i][j] = LU[i][j] * scale

    return LU, sp

# def get_splu(sparsity_pattern):
#     n = sparsity_pattern.shape[0]

#     zero_idxs = sparsity_pattern == 0

#     # Column ordering type
#     permc_spec = "MMD_AT_PLUS_A"
#     # permc_spec = "COLAMD"
#     # permc_spec = "NATURAL"

#     diag_pivot_thresh = 0.
#     options = {}
#     options["Equil"] = False # default True
#     options["RowPerm"] = "NOROWPERM"
#     options["PrintStat"] = True
#     options["SymmetricMode"] = True

#     M = np.random.rand(n, n)
#     M += 1
#     # Make diagonally dominant just so splu actually finishes without error and so perm_c is equal to perm_r
#     np.fill_diagonal(M, M.diagonal() + 10 + 100*np.random.rand(n))
#     M[zero_idxs] = 0
#     M_sparse = sparse.csc_array(M)

#     M_sparse_lu = splu(M_sparse, permc_spec=permc_spec, diag_pivot_thresh=diag_pivot_thresh, options=options)
#     assert(np.all(M_sparse_lu.perm_c == M_sparse_lu.perm_r))

#     return M_sparse_lu

# def lu_get_sparsity(A):
#     n = A.shape[0]
#     assert(n == A.shape[1])

#     sp = np.full(A.shape, False) # Set to True if entry can be nonzero; False if entry is guaranteed to be nonzero
#     sp[A != 0] = True

#     def is_nonzero(sp_ij, A_ij, A_ik, A_kj):
#         return (sp_ij != 0 or A_ij or (A_ik and A_kj))

#     for j in range(n):
#         # U
#         for i in range(j):
#             for k in range(i):
#                 sp[i][j] = is_nonzero(sp[i][j], sp[i][j], sp[i][k], sp[k][j])
#                 sp[i][j] = sp[i][j] or (sp[i][k] * sp[k][j])

#         # diag, L
#         for i in range(j,n):
#             for k in range(j):
#                 sp[i][j] = is_nonzero(sp[i][j], sp[i][j], sp[i][k], sp[k][j])
#                 sp[i][j] = sp[i][j] or (sp[i][k] * sp[k][j])

#         # Only allow diagonal pivots
#         if sp[j][j] == False:
#             raise ValueError("j = {}, sp[j][j] = {}".format(j, sp[j][j]))

#     return sp



np.random.seed(1)
A = np.array([[2, 5, 8, 7], [5, 2, 2, 8], [7, 5, 6, 6], [5, 4, 4, 8]], dtype=float)

sparsity_pattern = np.load("sparsity_pattern.npy")
zero_idxs = sparsity_pattern == 0

n = sparsity_pattern.shape[0]

# Get column ordering
permc_spec = "MMD_AT_PLUS_A"
# permc_spec = "COLAMD"
# permc_spec = "NATURAL"

diag_pivot_thresh = 0.
options = {}
options["Equil"] = False # default True
options["RowPerm"] = "NOROWPERM"
options["PrintStat"] = True
options["SymmetricMode"] = True

# tmp = sparsity_pattern
# tmp = np.random.rand(n, n)
# tmp += 1
# # Make diagonally dominant just so splu actually finishes without error and so perm_c is equal to perm_r
# np.fill_diagonal(tmp, tmp.diagonal() + 10 + 100*np.random.rand(n))
# tmp[zero_idxs] = 0
# M_sparse = sparse.csc_array(tmp)

# M_sparse_lu = splu(M_sparse, permc_spec=permc_spec, options=options)
M_sparse_lu = get_splu(sparsity_pattern)
L = M_sparse_lu.L
U = M_sparse_lu.U
assert(np.all(M_sparse_lu.perm_c == M_sparse_lu.perm_r))

print("M_sparse_lu.L.nnz + M_sparse_lu.U.nnz - n = {}".format(M_sparse_lu.L.nnz + M_sparse_lu.U.nnz - n)) # diagonal is double counted
# print("# nonzeros of dense triangular matrix = {}".format(n * (n + 1) // 2))

# Reconstruct M
# Construct sparse (inverse) permutation matrices
# Note: Pr is equal to np.eye(n)[invert_permutation(perm_r), :] and similarly for Pc
Pr = sparse.csc_array((np.ones(n), (M_sparse_lu.perm_r, np.arange(n))))
Pc = sparse.csc_array((np.ones(n), (np.arange(n), M_sparse_lu.perm_c)))
M = (Pr.T @ (M_sparse_lu.L @ M_sparse_lu.U) @ Pc.T).toarray()

# Permute columns (and rows)
M_p = M[:, invert_permutation(M_sparse_lu.perm_c)][invert_permutation(M_sparse_lu.perm_r), :]
# sm_p = Pr @ sm @ Pc

M_lu, M_sp = lu_and_sparsity(M_p)
M_l = np.tril(M_lu, -1)
np.fill_diagonal(M_l, 1.)
M_u = np.triu(M_lu)
assert_close(M_p, M_l @ M_u)

assert_close((M_lu != 0).astype(int), (M_sp != 0).astype(int))

print("np.count_nonzero(M_lu) = {}".format(np.count_nonzero(M_lu)))
print("np.count_nonzero(M_l) = {}".format(np.count_nonzero(M_l)))
print("np.count_nonzero(M_u) = {}".format(np.count_nonzero(M_u)))

M_sp2 = get_lu_sparsity(M_p)
assert_close((M_lu != 0).astype(int), (M_sp2 != 0).astype(int))

sp_permuted = get_lu_sparsity(sparsity_pattern[:, invert_permutation(M_sparse_lu.perm_c)][invert_permutation(M_sparse_lu.perm_r), :])
assert_close((sp_permuted != 0).astype(int), (M_sp2 != 0).astype(int))

# #
# # A = sparse.csc_array(sparsity_pattern, dtype=float)

# A = np.random.rand(sparsity_pattern.shape[0], sparsity_pattern.shape[1])
# A += 1
# low = 1
# high = 10
# # A = np.random.randint(low, high, size=sparsity_pattern.shape).astype(float) # integers between low and high-1

# # Make more diagonally dominant
# np.fill_diagonal(A, A.diagonal() + high)

# A[zero_idxs] = 0.


p, l, u = lu(A)
np.testing.assert_allclose(A, p @ l @ u)
np.testing.assert_allclose(p.T @ A, l @ u)

p_idx, _, _ = lu(A, p_indices=True)
p_idx_inv = invert_permutation(p_idx)
np.testing.assert_allclose(A, l[p_idx, :] @ u)
np.testing.assert_allclose(A[p_idx_inv, :], l @ u)


l1 = l; u1 = u
# A_LU, _ = lu_and_sparsity(A[p_idx_inv, :])
A_LU, _ = lu_and_sparsity(A)
l2 = np.tril(A_LU, -1)
np.fill_diagonal(l2, 1)
u2 = np.triu(A_LU)
assert_close(A, l2 @ u2)


breakpoint()