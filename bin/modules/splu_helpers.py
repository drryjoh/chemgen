import numpy as np
from scipy.linalg import lu
from scipy import sparse
from scipy.sparse.linalg import splu, spilu

def invert_permutation(p):
    """Return an array s with which np.array_equal(arr[p][s], arr) is True.
    The array_like argument p must be some permutation of 0, 1, ..., len(p)-1.

    # https://stackoverflow.com/questions/11649577/how-to-invert-a-permutation-array-in-numpy
    """
    p = np.asanyarray(p) # in case p is a tuple, etc.
    s = np.empty_like(p)
    s[p] = np.arange(p.size)
    return s

def make_diagonally_dominant(M):
    np.fill_diagonal(M, M.diagonal() + 10 + 100*np.random.rand(M.shape[0]))

    return M

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
    M = make_diagonally_dominant(M)
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

def get_splu_sparsity_empirical(A):
    n = A.shape[0]
    assert(n == A.shape[1])

    zero_idxs = A == 0
    nonzero_idxs = A != 0

    # Column ordering type
    # permc_spec = "MMD_AT_PLUS_A"
    # permc_spec = "COLAMD"
    permc_spec = "NATURAL"

    diag_pivot_thresh = 0.
    options = {}
    options["Equil"] = False # default True
    options["RowPerm"] = "NOROWPERM"
    options["PrintStat"] = True
    options["SymmetricMode"] = True

    n_test = 10 # number of test matrices
    for i in range(n_test):
        print("Test matrix {}".format(i+1))
        M = np.random.rand(n, n)
        M += 1
        # Make diagonally dominant just so splu actually finishes without error and so perm_c is equal to perm_r
        M = make_diagonally_dominant(M)
        M[zero_idxs] = 0
        M_sparse = sparse.csc_array(M)

        M_sparse_lu = splu(M_sparse, permc_spec=permc_spec, diag_pivot_thresh=diag_pivot_thresh, options=options)
        assert(np.all(M_sparse_lu.perm_c == M_sparse_lu.perm_r))

        L = M_sparse_lu.L
        U = M_sparse_lu.U

        if i == 0:
            # store indices
            L_indices = L.indices
            L_indptr = L.indptr
            U_indices = U.indices
            U_indptr = U.indptr
        else:
            # compare
            np.testing.assert_array_equal(L_indices, L.indices)
            np.testing.assert_array_equal(L_indptr, L.indptr)
            np.testing.assert_array_equal(U_indices, U.indices)
            np.testing.assert_array_equal(U_indptr, U.indptr)

    # Convert to sparsity pattern (stored as dense matrix of type int)
    L.data[:] = 1
    L.setdiag(0)
    U.data[:] = 1
    sp = (L.toarray() + U.toarray()).astype(int)
    
    return sp