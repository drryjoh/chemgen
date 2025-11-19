#!python3
import numpy as np, os, glob
d1 = 1/np.sort(np.load( "data/eigs/eigs_cgc_36.npy"))
d2 = 1/np.sort(np.load( "data/eigs/eigs_cgt_36.npy"))
for d1i, d2i in zip(d1, d2):
    print(d1i, d2i)


