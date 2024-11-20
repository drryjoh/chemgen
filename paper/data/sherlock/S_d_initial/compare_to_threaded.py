import numpy as np
import matplotlib.pyplot as plt

#load threaded info
#H2
t_h2_tbb = np.load("../t1_initial/data/chemgen_h2_16_best_time.npy")
pts_h2_tbb = np.load("../t1_initial/data/chemgen_h2_16_pts.npy")
times_h2_mpi = np.load("./data/chemgen_h2_times.npy")
ranks_h2_mpi = np.load("./data/chemgen_h2_ranks.npy")
pts_h2_mpi = np.load("./data/chemgen_h2_pts.npy")

t_ffcm2_tbb = np.load("../t1_initial/data/chemgen_ffcm2_16_best_time.npy")
pts_ffcm2_tbb = np.load("../t1_initial/data/chemgen_ffcm2_16_pts.npy")
times_ffcm2_mpi = np.load("./data/chemgen_ffcm2_times.npy")
ranks_ffcm2_mpi = np.load("./data/chemgen_ffcm2_ranks.npy")
pts_ffcm2_mpi = np.load("./data/chemgen_ffcm2_pts.npy")

ranks = np.load("./data/chemgen_ffcm2_ranks.npy")
for i, r in enumerate(ranks[0,:]):
    plt.figure()
    #plt.loglog(pts_h2_mpi, times_h2_mpi[:,i],'ok', label = "S-d")
    #plt.loglog(pts_h2_tbb, t_h2_tbb,'ok',label = 'T1')
    plt.semilogx(pts_h2_tbb, times_h2_mpi[:,i]/t_h2_tbb,'ok',label = 'FFCM2 H2 model')

    #plt.loglog(pts_ffcm2_mpi, times_ffcm2_mpi[:,i],'ok', label = "S-d")
    #plt.loglog(pts_ffcm2_tbb, t_ffcm2_tbb,'or',label = 'T1')
    plt.semilogx(pts_ffcm2_tbb, times_ffcm2_mpi[:,i]/t_ffcm2_tbb,'or',label = 'FFCM2 Full model')
    plt.semilogx(pts_ffcm2_tbb, t_ffcm2_tbb*0+1.0,'--k',label = '')

    plt.legend()
    plt.ylim([0,20])
    plt.xlabel("Number of points")
    plt.ylabel(f"Compute time: resource $n_{{proc}} = {int(r)}$ $n_{{threads}}={int(16/r)}$")
plt.show()
