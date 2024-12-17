import numpy as np
import matplotlib.pyplot as plt
ranks =[8]# [1,2,4,8]
file_names = [f"chemgen_ffcm2_{x}.out" for x in ranks]

color_idx = np.linspace(0,1,len(ranks))
colors = []
[colors.append(plt.cm.jet(idx)) for idx in color_idx]


for x in ranks:
    file_names.append(f"../t2_initial/chemgen_ffcm2_{x}.out")
[colors.append(plt.cm.jet(idx)) for idx in color_idx]

for x in ranks:
    file_names.append(f"../t3_initial/chemgen_ffcm2_{x}.out")
[colors.append(plt.cm.jet(idx)) for idx in color_idx]

labels = []

mpi_16_pts = np.load("../S_d_initial/data/chemgen_ffcm2_pts.npy")
mpi_16_times = np.load("../S_d_initial/data/chemgen_ffcm2_times.npy")
thread_16_times  = np.load("../t1_initial/data/chemgen_ffcm2_best_time.npy")
baseline = mpi_16_times[:,-1]
for k, file in enumerate(file_names):
    pts = []
    times = []
    n_threads = []
    marker = 'o'
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
            if "t2_initial" in file:
                labels.append(f"T2: Threads = {threads}\n Ranks = {int(16/threads)} ")
                marker = 'o'
            elif "t3_initial" in file:
                labels.append(f"T3a: Threads = {threads}\n Ranks = {int(16/threads)} ")
                marker = 'o'
            else: 
                labels.append(f"T3b: Threads = {threads}\n Ranks = {int(16/threads)} ")
                marker = 'x'
        
        if len(pts) == len(times): 
            plt.semilogx(pts, baseline/times,f'-{marker}',mfc="white",label = labels[k], color = colors[k])
        else:
            plt.semilogx(pts[0:len(times)], baseline[0:len(times)]/times,'-o',mfc="white",label = labels[k], color = colors[k])
plt.legend(fontsize = 8, ncol=2, markerscale=.75)
plt.ylabel("Speed up over pure MPI on 16 cores")
plt.xlabel("Number of Degrees of Freedom")
plt.savefig("speedup.png", dpi=300)

plt.show()
