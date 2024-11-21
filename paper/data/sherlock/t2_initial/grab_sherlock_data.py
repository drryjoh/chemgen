import numpy as np
import matplotlib.pyplot as plt
file_names = ["chemgen_ffcm2_2.out", "chemgen_ffcm2_4.out", "chemgen_ffcm2_8.out", "chemgen_ffcm2_1.out"]#, "chemgen_h2.out", "chemgen_h2_16.out", "chemgen_ffcm2_16.out"]
file_names = [f"chemgen_ffcm2_{x}.out" for x in [1,2,4,8,16]]
labels = []

mpi_16_pts = np.load("../S_d_initial/data/chemgen_ffcm2_pts.npy")
mpi_16_times = np.load("../S_d_initial/data/chemgen_ffcm2_times.npy")
thread_16_times  = np.load("../t1_initial/data/chemgen_ffcm2_best_time.npy")
baseline = mpi_16_times[:,-1]
for k, file in enumerate(file_names):
    pts = []
    times = []
    n_threads = []
    with open(file,"r") as f:
        for line in f:
            if "Number of points" in line:
                pts.append(float(line.strip('\n').split('points')[1].replace(' ','')))
            if "time" in line:
                times.append(float(line.split(':')[1].split("seconds")[0].replace(' ','')))
            if "threads" in line:
                n_threads.append( line.split(':')[1].replace(' ','').replace('\n',''))
        if all(t == n_threads[0] for t in n_threads):
            threads = int(n_threads[0])
            labels.append(f"T2: Threads = {threads}\n Ranks = {int(16/threads)} ")

        if len(pts) == len(times): 
            plt.semilogx(pts, baseline/times,'-o',mfc="white",label = labels[k])
        else:
            plt.semilogx(pts[0:len(times)], baseline[0:len(times)]/times,'-o',mfc="white",label = labels[k])
plt.semilogx(mpi_16_pts, baseline /thread_16_times, 'ok', label="T1: Threads = 16")
plt.legend(fontsize = 6, ncol=2, markerscale=.75)
plt.ylabel("Speed up over pure MPI on 16 cores")
plt.show()
