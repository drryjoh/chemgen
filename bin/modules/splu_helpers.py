import numpy as np
from scipy.linalg import lu
from scipy import sparse
from scipy.sparse.linalg import splu, spilu

def get_splu(sparsity_pattern):
    n = sparsity_pattern.shape[0]

    zero_idxs = sparsity_pattern == 0

    # Column ordering type
    permc_spec = "MMD_AT_PLUS_A"
    # permc_spec = "COLAMD"
    # permc_spec = "NATURAL"

    diag_pivot_thresh = 0.
    options = {}
    options["Equil"] = False # default True
    options["RowPerm"] = "NOROWPERM"
    options["PrintStat"] = True
    options["SymmetricMode"] = True

    M = np.random.rand(n, n)
    M += 1
    # Make diagonally dominant just so splu actually finishes without error and so perm_c is equal to perm_r
    np.fill_diagonal(M, M.diagonal() + 10 + 100*np.random.rand(n))
    M[zero_idxs] = 0
    M_sparse = sparse.csc_array(M)

    M_sparse_lu = splu(M_sparse, permc_spec=permc_spec, diag_pivot_thresh=diag_pivot_thresh, options=options)
    assert(np.all(M_sparse_lu.perm_c == M_sparse_lu.perm_r))

    return M_sparse_lu

def get_lu_sparsity(A):
    n = A.shape[0]
    assert(n == A.shape[1])

    sp = np.full(A.shape, False) # Set to True if entry can be nonzero; False if entry is guaranteed to be nonzero
    sp[A != 0] = True

    def is_nonzero(sp_ij, A_ij, A_ik, A_kj):
        return (sp_ij != 0 or A_ij or (A_ik and A_kj))

    for j in range(n):
        # U
        for i in range(j):
            for k in range(i):
                sp[i][j] = is_nonzero(sp[i][j], sp[i][j], sp[i][k], sp[k][j])
                sp[i][j] = sp[i][j] or (sp[i][k] * sp[k][j])

        # diag, L
        for i in range(j,n):
            for k in range(j):
                sp[i][j] = is_nonzero(sp[i][j], sp[i][j], sp[i][k], sp[k][j])
                sp[i][j] = sp[i][j] or (sp[i][k] * sp[k][j])

        # Only allow diagonal pivots
        if sp[j][j] == False:
            raise ValueError("j = {}, sp[j][j] = {}".format(j, sp[j][j]))

    return sp