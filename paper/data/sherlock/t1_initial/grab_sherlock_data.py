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
with open("chemgen.out","r") as f:
    for line in f:
        if "Number of points" in line:
            pts.append(float(line.strip('\n').split('points')[1].replace(' ','')))
        if "point_reactions" in line:
            time = float(line.strip('\n').split(':')[1].split('seconds')[0])
            if 'Serial' in line:
                forward_reactions_serial.append(time)
            elif '20' in line:
                forward_reactions_20.append(time)
            elif '100' in line:
                print(time)
                forward_reactions_100.append(time)
            else:
                forward_reactions_no_chunk.append(time)
        if "point_progress_rates" in line:
            time = float(line.strip('\n').split(':')[1].split('seconds')[0])
            if 'Serial' in line:
                point_progress_rates_serial.append(time)
            elif '20' in line:
                point_progress_rates_20.append(time)
            elif '100' in line:
                print(time)
                point_progress_rates_100.append(time)
            else:
                point_progress_rates_no_chunk.append(time)
pts = np.array(pts)
number_of_passed = len(forward_reactions_serial)
pts = pts[:number_of_passed]
forward_reactions_serial = np.array(forward_reactions_serial)
forward_reactions_no_chunk = np.array(forward_reactions_no_chunk) 
forward_reactions_20 = np.array(forward_reactions_20) 
forward_reactions_100 = np.array(forward_reactions_100)

point_progress_rates_serial = np.array(point_progress_rates_serial)
point_progress_rates_no_chunk = np.array(point_progress_rates_no_chunk) 
point_progress_rates_20 = np.array(point_progress_rates_20) 
point_progress_rates_100 = np.array(point_progress_rates_100)

plt.plot(pts, forward_reactions_serial/forward_reactions_no_chunk, 'ok')
plt.plot(pts[:-1], point_progress_rates_serial[:-1]/point_progress_rates_100, 'or')
plt.show()
