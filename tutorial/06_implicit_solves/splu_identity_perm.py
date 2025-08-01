import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import splu, spilu

# Load sparsity pattern
sp_perm = np.load("input_sparsity.npy")
n = sp_perm.shape[0]

''' BEGIN: REPLACE A BELOW WITH ACTUAL MATRIX '''
np.random.seed(1)
A = np.random.rand(n, n) # dense matrix format
A[sp_perm == 0] = 0
''' END '''

# Create triplets
row, col = np.nonzero(sp_perm)
data = A[row, col]

# Create sparse matrix
A_sparse = sparse.csc_array((data, (row, col)), shape=(n, n), dtype=float)

# Perform sparse LU without any permutations
options = {}
options["Equil"] = False # default True
options["RowPerm"] = "NOROWPERM"
options["SymmetricMode"] = False

lu = splu(A_sparse, permc_spec="NATURAL", diag_pivot_thresh=0., options=options)

# Check that splu permutation is identity
assert(np.all(lu.perm_c == lu.perm_r))
assert(np.all(lu.perm_c == np.arange(n)))

# Check that lu reproduces A
assert(np.linalg.norm(A - lu.L @ lu.U) < 1e-10)
