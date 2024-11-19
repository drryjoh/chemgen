import numpy as np
import matplotlib.pyplot as plt
pts = []
forward_reactions_serial = []
forward_reactions_no_chunk = []
forward_reactions_20 = []
forward_reactions_100 = []

point_progress_rates_serial = []
point_progress_rates_no_chunk = []
point_progress_rates_20 = []
point_progress_rates_100 = []

source_serial = []
source_no_chunk = []
source_20 = []
source_100 = []
with open("chemgen_h2.out","r") as f:
    for line in f:
        if "Number of points" in line:
            pts.append(float(line.strip('\n').split('points')[1].replace(' ','')))
        if "point_reactions" in line:
            time = float(line.strip('\n').split(':')[1].split('seconds')[0])
            if 'Serial' in line:
                forward_reactions_serial.append(time)
            elif 'chunk 20:' in line:
                forward_reactions_20.append(time)
            elif 'chunk 100:' in line:
                forward_reactions_100.append(time)
            else:
                forward_reactions_no_chunk.append(time)
        if "point_progress_rates" in line:
            time = float(line.strip('\n').split(':')[1].split('seconds')[0])
            if 'Serial' in line:
                point_progress_rates_serial.append(time)
            elif 'chunk 20:' in line:
                point_progress_rates_20.append(time)
            elif 'chunk 100:' in line:
                point_progress_rates_100.append(time)
            else:
                point_progress_rates_no_chunk.append(time)
        if "species_source" in line:
            time = float(line.strip('\n').split(':')[1].split('seconds')[0])
            if 'Serial' in line:
                source_serial.append(time)
            elif 'chunk 20:' in line:
                source_20.append(time)
            elif 'chunk 100:' in line:
                source_100.append(time)
            else:
                source_no_chunk.append(time)
pts = np.array(pts)
number_of_passed = len(forward_reactions_serial)
pts = pts[:number_of_passed]
forward_reactions_serial = np.array(forward_reactions_serial)
forward_reactions_no_chunk = np.array(forward_reactions_no_chunk) 
forward_reactions_20 = np.array(forward_reactions_20) 
forward_reactions_100 = np.array(forward_reactions_100)

tbb_fr = np.array([forward_reactions_no_chunk, forward_reactions_20, forward_reactions_100])
best_time_ft  = np.min(tbb_fr, axis=0)

point_progress_rates_serial = np.array(point_progress_rates_serial)
point_progress_rates_no_chunk = np.array(point_progress_rates_no_chunk) 
point_progress_rates_20 = np.array(point_progress_rates_20) 
point_progress_rates_100 = np.array(point_progress_rates_100)

tbb_progress = np.array([point_progress_rates_no_chunk, point_progress_rates_20, point_progress_rates_100])
best_time_progress  = np.min(tbb_progress, axis=0)

source_serial = np.array(source_serial)
source_no_chunk = np.array(source_no_chunk) 
source_20 = np.array(source_20) 
source_100 = np.array(source_100)


tbb_source = np.array([source_no_chunk, source_20, source_100])
best_time_source  = np.min(tbb_source, axis=0)

best_time = best_time_ft + best_time_progress + best_time_source
serial_time = forward_reactions_serial + point_progress_rates_serial + source_serial

plt.figure()
plt.semilogx(pts, forward_reactions_serial/forward_reactions_no_chunk, '-ok',label = "$k_f$ Default TBB")
plt.semilogx(pts, forward_reactions_serial/forward_reactions_20, '-ob',label = "$k_f$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/20$")
plt.semilogx(pts, forward_reactions_serial/forward_reactions_100, '-or',label = "$k_f$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
plt.legend()
plt.xlabel("# of points")
plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")


plt.figure()
plt.semilogx(pts, point_progress_rates_serial/point_progress_rates_no_chunk, '-^k',label = "$q_i$ Default TBB")
plt.semilogx(pts, point_progress_rates_serial/point_progress_rates_20, '-^b',label = "$q_i$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/20$")
plt.semilogx(pts, point_progress_rates_serial/point_progress_rates_100, '-^r',label = "$q_i$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
plt.legend()
plt.xlabel("# of points")
plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")


plt.figure()
plt.semilogx(pts, source_serial/source_no_chunk, '-xk', label="$\\omega_i$ default TBB")
plt.semilogx(pts, source_serial/source_20, '-xb', label="$\\omega_i$ \n $n_{chunk} = n_{r}\\times n_{pts}/20$")
plt.semilogx(pts, source_serial/source_100, '-xr', label="$\\omega_i$ \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
plt.legend()
plt.xlabel("# of points")
plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")



plt.figure()

plt.semilogx(pts, serial_time/best_time, '-k', label="$\\omega_i$ \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
plt.show()
