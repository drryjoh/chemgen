import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import splu, spilu
import scikits.umfpack as um

def invert_permutation(p):
    """Return an array s with which np.array_equal(arr[p][s], arr) is True.
    The array_like argument p must be some permutation of 0, 1, ..., len(p)-1.

    # https://stackoverflow.com/questions/11649577/how-to-invert-a-permutation-array-in-numpy
    """
    p = np.asanyarray(p) # in case p is a tuple, etc.
    s = np.empty_like(p)
    s[p] = np.arange(p.size)
    return s

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

np.random.seed(1)
# m = 10
# p = 0.01
# A = sparse.random(m,m,p)
# f = splu(A, options={"ColPerm":"NATURAL"}) # no difference
# print(assert_close(f.perm_c,np.arange(m)))
# print(f.perm_c[:10])

eps = np.finfo(float).eps

# n = 10
# A = sparse.identity(n, format="csc")
# f = splu(A, permc_spec="NATURAL")
# print(A)

sparsity_pattern = np.load("sparsity_pattern.npy")
zero_idxs = sparsity_pattern == 0

#
# A = sparse.csc_array(sparsity_pattern, dtype=float)

A = np.random.rand(sparsity_pattern.shape[0], sparsity_pattern.shape[1])
A += 1
low = 1
high = 100
# A = np.random.randint(low, high, size=sparsity_pattern.shape).astype(float) # integers between low and high-1

n = A.shape[0]

# Make more diagonally dominant
np.fill_diagonal(A, A.diagonal() + np.arange(high, high+n))

A[zero_idxs] = 0.

n_nonzeros = np.count_nonzero(sparsity_pattern)

row = np.zeros(n_nonzeros, dtype=np.int32)
col = np.zeros(n_nonzeros, dtype=np.int32)
data = np.zeros(n_nonzeros)

k = 0
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        if sparsity_pattern[i][j] != 0:
            row[k] = i
            col[k] = j
            data[k] = A[i][j]
            k += 1

np.testing.assert_equal(k, n_nonzeros)

# sm = sparse.csc_array(A)
sm = sparse.csc_array((data, (row, col)), shape=(n, n), dtype=float)

# dm = np.ones([5,5])
# # dm[0,0] = 1.
# dm[0,0] = eps
# dm[1,1] = eps
# sm = sparse.csc_array(dm)

# triplets for sparse matrix
# row = np.array([0, 2, 2, 0, 1, 2, 0])
# col = np.array([0, 0, 1, 2, 2, 2, 1])
# data = np.array([1, 2, 3, 4, 5, 6, 1])

# end = -1

# # row = row[:end]; col = col[:end]; data = data[:end]

# n = 3
# sm = sparse.csc_array((data, (row, col)), shape=(n, n), dtype=float)

permc_spec = "MMD_AT_PLUS_A"
# permc_spec = "COLAMD"
# permc_spec = "NATURAL"

diag_pivot_thresh = 0.
options = {}
options["Equil"] = False # default True
options["RowPerm"] = "NOROWPERM"
options["PrintStat"] = True
options["SymmetricMode"] = False
# options["ParSymbFact"] = "YES"

lu = splu(sm, permc_spec=permc_spec, diag_pivot_thresh=diag_pivot_thresh, options=options)
# lu = spilu(sm, permc_spec=permc_spec, diag_pivot_thresh=diag_pivot_thresh, options=options, drop_tol=0., fill_factor=100, drop_rule="basic")

tmp = sparsity_pattern
tmp[:,:] = 1
np.fill_diagonal(tmp, tmp.diagonal() + 100)
tmp[zero_idxs] = 0
sm_sparsity = sparse.csc_array(tmp)
lu_sparsity = splu(sm_sparsity, permc_spec=permc_spec, diag_pivot_thresh=diag_pivot_thresh, options=options)
assert_close(lu_sparsity.perm_c, lu.perm_c)
# assert_close(lu_sparsity.perm_r, lu.perm_r) # this can be different
# print("L = \n{}".format(lu.L))
# print("U = \n{}".format(lu.U))

print("perm_c = {}".format(lu.perm_c))
print("perm_r = {}".format(lu.perm_r))

# print("np.amin(lu.perm_c-lu.perm_r) = {}".format(np.amin(lu.perm_c-lu.perm_r)))

assert(np.all(lu.perm_c == lu.perm_r))

# Construct sparse (inverse) permutation matrices
# Note: Pr is equal to np.eye(n)[invert_permutation(lu.perm_r), :] and similarly for Pc
Pr = sparse.csc_array((np.ones(n), (lu.perm_r, np.arange(n))))
Pc = sparse.csc_array((np.ones(n), (np.arange(n), lu.perm_c)))
m_check = (Pr.T @ (lu.L @ lu.U) @ Pc.T).toarray()
assert_close(A, m_check)
assert_close(Pr @ sm @ Pc, lu.L @ lu.U)
# print("m_check = \n{}".format(m_check))


# Recompute with perm_c and perm_r already included before splu call
A_p = A[:, invert_permutation(lu.perm_c)][invert_permutation(lu.perm_r), :]
sm_p = Pr @ sm @ Pc
# print("np.linalg.norm(A_p - sm_p) = {}".format(np.linalg.norm(A_p - sm_p)))
lu_p = splu(sm_p, permc_spec="NATURAL", diag_pivot_thresh=diag_pivot_thresh, options=options)
assert(lu.L.size == lu_p.L.size)
assert(lu.U.size == lu_p.U.size)
assert(np.all(lu_p.perm_c == lu_p.perm_r))
assert(np.all(lu_p.perm_c == np.arange(n)))
assert_close(A_p, lu_p.L @ lu_p.U)
# print("np.linalg.norm(A_p - (lu_p.L @ lu_p.U).toarray()) = {}".format(np.linalg.norm(A_p - (lu_p.L @ lu_p.U).toarray())))
# print("np.linalg.norm((sm_p - lu_p.L @ lu_p.U).toarray()) = {}".format(np.linalg.norm((sm_p - lu_p.L @ lu_p.U).toarray())))

# Try dense lu
P, L, U = scipy.linalg.lu(A_p)
assert(np.amin(P.diagonal()) == 1.)
assert_close(A_p, L @ U)
# print("np.linalg.norm(A_p - L @ U) = {}".format(np.linalg.norm(A_p - L @ U)))
assert(lu.L.size == np.count_nonzero(L))
assert(lu.U.size == np.count_nonzero(U))

# U_arr = lu.U.toarray()
# Pc_arr = Pc.toarray()
# U_times_PcT_arr = (lu.U @ Pc.T).toarray()

# number of zeros in sparse matrix
n_zeros_sparse = sm.size - sm.count_nonzero()

print("n_zeros_sparse = {}".format(n_zeros_sparse))
print("sm.count_nonzero() = {}".format(sm.count_nonzero()))
print("sm.size = {}".format(sm.size))

# np.count_nonzero(np.tril(np.ones([n,n])))

print("# nonzeros of dense triangular matrix = {}".format(n * (n + 1) // 2))
print("lu.L.size = {}".format(lu.L.size))
print("lu.U.size = {}".format(lu.U.size))
print("lu.L.nnz = {}".format(lu.L.nnz))
print("lu.U.nnz = {}".format(lu.U.nnz))

print("lu.L.count_nonzero() = {}".format(lu.L.count_nonzero()))
print("lu.U.count_nonzero() = {}".format(lu.U.count_nonzero()))


Adiff = A - m_check

print("np.linalg.norm(Adiff) = {}".format(np.linalg.norm(Adiff)))


# https://scicomp.stackexchange.com/questions/3229/permute-a-matrix-in-place-in-numpy
# https://www.reddit.com/r/learnpython/comments/sg8w28/automatically_generate_square_permutation_matrices/

# Contruct the solver.
umf = um.UmfpackContext() # Use default 'di' family of UMFPACK routines.
umf.control[um.UMFPACK_PRL] = 4 # Let's be more verbose.
# mtx0 = sparse.csc_array((data, (row, col)), shape=(n, n), dtype=np.dtype(np.int32))
# breakpoint()
umf.symbolic(sm)
# umf.report_symbolic()

umf.funs.save_symbolic(umf.symbolic(sm), "/Users/eching/Codes/chemgen/tutorial/06_implicit_solves/symbolic.txt")

breakpoint()