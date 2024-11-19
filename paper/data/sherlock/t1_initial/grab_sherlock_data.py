import numpy as np
import matplotlib.pyplot as plt

file_names = ["chemgen_ffcm2.out", "chemgen_h2.out"]#,"chemgen_h2.out"]
for file in file_names:
    pts = []
    forward_reactions_serial = []
    forward_reactions_no_chunk = []
    forward_reactions_10nthreads = []
    forward_reactions_nthreads = []

    point_progress_rates_serial = []
    point_progress_rates_no_chunk = []
    point_progress_rates_10nthreads = []
    point_progress_rates_nthreads = []

    source_serial = []
    source_no_chunk = []
    source_10nthreads = []
    source_nthreads = []
    prepend = "data/{file}".format(file = file.split(',')[0])
    
    with open(file,"r") as f:
        for line in f:
            if "Number of points" in line:
                pts.append(float(line.strip('\n').split('points')[1].replace(' ','')))
            if "point_reactions" in line:
                time = float(line.strip('\n').split(':')[1].split('seconds')[0])
                if 'Serial' in line:
                    forward_reactions_serial.append(time)
                elif 'chunk 10*n_threads:' in line:
                    forward_reactions_10nthreads.append(time)
                elif 'chunk n_threads:' in line:
                    forward_reactions_nthreads.append(time)
                else:
                    forward_reactions_no_chunk.append(time)
            if "point_progress_rates" in line:
                time = float(line.strip('\n').split(':')[1].split('seconds')[0])
                if 'Serial' in line:
                    point_progress_rates_serial.append(time)
                elif 'chunk 10*n_threads:' in line:
                    point_progress_rates_10nthreads.append(time)
                elif 'chunk n_threads:' in line:
                    point_progress_rates_nthreads.append(time)
                else:
                    point_progress_rates_no_chunk.append(time)
            if "species_source" in line:
                time = float(line.strip('\n').split(':')[1].split('seconds')[0])
                if 'Serial' in line:
                    source_serial.append(time)
                elif 'chunk 10*n_threads:' in line:
                    source_10nthreads.append(time)
                elif 'chunk n_threads:' in line:
                    source_nthreads.append(time)
                else:
                    source_no_chunk.append(time)
    pts = np.array(pts)
    np.save(f"{prepend}_pts.npy", pts)

    forward_reactions_serial = np.array(forward_reactions_serial)
    forward_reactions_no_chunk = np.array(forward_reactions_no_chunk) 
    forward_reactions_10nthreads = np.array(forward_reactions_10nthreads) 
    forward_reactions_nthreads = np.array(forward_reactions_nthreads)
    np.save(f"{prepend}_forward_reactions_serial.npy", forward_reactions_serial)
    np.save(f"{prepend}_forward_reactions_no_chunk.npy", forward_reactions_no_chunk)
    np.save(f"{prepend}_forward_reactions_10nthreads.npy", forward_reactions_10nthreads)
    np.save(f"{prepend}_forward_reactions_nthreads.npy", forward_reactions_nthreads)

    tbb_fr = np.array([forward_reactions_no_chunk, forward_reactions_10nthreads, forward_reactions_nthreads])
    best_time_ft  = np.min(tbb_fr, axis=0)
    np.save(f"{prepend}_best_time_ft.npy", pts)

    point_progress_rates_serial = np.array(point_progress_rates_serial)
    point_progress_rates_no_chunk = np.array(point_progress_rates_no_chunk) 
    point_progress_rates_10nthreads = np.array(point_progress_rates_10nthreads) 
    point_progress_rates_nthreads = np.array(point_progress_rates_nthreads)
    np.save(f"{prepend}_point_progress_rates_serial.npy", point_progress_rates_serial)
    np.save(f"{prepend}_point_progress_rates_no_chunk.npy", point_progress_rates_no_chunk)
    np.save(f"{prepend}_point_progress_rates_10nthreads.npy", point_progress_rates_10nthreads)
    np.save(f"{prepend}_point_progress_rates_nthreads.npy", point_progress_rates_nthreads)


    tbb_progress = np.array([point_progress_rates_no_chunk, point_progress_rates_10nthreads, point_progress_rates_nthreads])
    best_time_progress  = np.min(tbb_progress, axis=0)
    np.save(f"{prepend}_best_time_progress.npy", best_time_progress)

    source_serial = np.array(source_serial)
    source_no_chunk = np.array(source_no_chunk) 
    source_10nthreads = np.array(source_10nthreads) 
    source_nthreads = np.array(source_nthreads)
    np.save(f"{prepend}_source_serial.npy", source_serial)
    np.save(f"{prepend}_source_no_chunk.npy", source_no_chunk)
    np.save(f"{prepend}_source_10nthreads.npy", source_10nthreads)
    np.save(f"{prepend}_source_nthreads.npy", source_nthreads)

    tbb_source = np.array([source_no_chunk, source_10nthreads, source_nthreads])
    best_time_source  = np.min(tbb_source, axis=0)
    np.save(f"{prepend}_best_time_source.npy", best_time_source)

    best_time = best_time_ft + best_time_progress + best_time_source
    serial_time = forward_reactions_serial + point_progress_rates_serial + source_serial

    plt.figure()
    plt.semilogx(pts, forward_reactions_serial/forward_reactions_no_chunk, '-ok',label = "$k_f$ Default TBB")
    plt.semilogx(pts, forward_reactions_serial/forward_reactions_10nthreads, '-ob',label = "$k_f$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/(10n_{threads})$")
    plt.semilogx(pts, forward_reactions_serial/forward_reactions_nthreads, '-or',label = "$k_f$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
    plt.legend()
    plt.xlabel("# of points")
    plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")


    plt.figure()
    plt.semilogx(pts, point_progress_rates_serial/point_progress_rates_no_chunk, '-^k',label = "$q_i$ Default TBB")
    plt.semilogx(pts, point_progress_rates_serial/point_progress_rates_10nthreads, '-^b',label = "$q_i$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/(10n_{threads})$")
    plt.semilogx(pts, point_progress_rates_serial/point_progress_rates_nthreads, '-^r',label = "$q_i$ TBB \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
    plt.legend()
    plt.xlabel("# of points")
    plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")


    plt.figure()
    plt.semilogx(pts, source_serial/source_no_chunk, '-xk', label="$\\omega_i$ default TBB")
    plt.semilogx(pts, source_serial/source_10nthreads, '-xb', label="$\\omega_i$ \n $n_{chunk} = n_{r}\\times n_{pts}/(10n_{threads})$")
    plt.semilogx(pts, source_serial/source_nthreads, '-xr', label="$\\omega_i$ \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
    plt.legend()
    plt.xlabel("# of points")
    plt.ylabel("Speed up over serial, $\\tau_{S}/\\tau_{T1}$")



    plt.figure()
    plt.semilogx(pts, serial_time/best_time, '-k', label="$\\omega_i$ \n $n_{chunk} = n_{r}\\times n_{pts}/n_{threads}$")
plt.show()
