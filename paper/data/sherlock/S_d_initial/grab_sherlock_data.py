import numpy as np
import matplotlib.pyplot as plt
file_names  = ["chemgen_h2.out", "chemgen_ffcm2.out"]

for file in file_names:
    pts = []
    ranks_per_points = []
    times_per_points = []
    ranks = []
    times = []
    current_rank  = 0
    time = 0
    prepend = "data/{file}".format(file = file.split('.')[0])
    
    with open(file,"r") as f:
        for line in f:
            if "Number of points" in line:
                pts.append(float(line.rstrip('\n').split('points')[1].replace(' ','')))
            if "MPI processes" in line:
                rank = float(line.split('MPI processes')[0].split("with")[1].replace(' ',''))
                if rank < current_rank:
                    # we restarted
                    ranks_per_points.append(ranks)
                    current_rank = rank
                    ranks = []
                    ranks.append(rank)
                else:
                    ranks.append(rank)
                    current_rank = rank
            if "Rank 1" in line:
                time  = float(line.split(':')[1].split('seconds')[0].replace(' ',''))
                times.append(time)
    ranks_per_points.append(ranks)
    ranks_per_points = np.array(ranks_per_points)
    shape = np.shape(ranks_per_points)
    times_per_points = np.zeros(shape)
    for i, ranks in enumerate(ranks_per_points):
        for j, rank in enumerate(ranks):
            times_per_points[i,j] = times[j+shape[1]*i]
    plt.figure()
    color_idx = np.linspace(0,1,shape[0])
    for k, (ranks, times) in enumerate(zip(ranks_per_points, times_per_points)):
        plt.plot(ranks, times[0]/times, label = f"$n_{{p}} = {int(pts[k])}$", color=plt.cm.jet(color_idx[k]))
    plt.plot(ranks, ranks,'--k', label = "ideal scaling") 
    plt.legend()
    plt.ylim([0,20])
    plt.xlabel("Ranks")
    plt.ylabel("Speed up $\\tau_{S}/\\tau_{S-d}$")
    
    np.save(f"{prepend}_times.npy", times_per_points)
    np.save(f"{prepend}_ranks.npy", ranks_per_points)
    np.save(f"{prepend}_pts.npy", pts)
plt.show()
