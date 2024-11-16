import numpy as np
import matplotlib.pyplot as plt
pts = []
forward_reactions_serial = []
forward_reactions_no_chunk = []
forward_reactions_20 = []
forward_reactions_100 = []
with open("test_job.55647875.out","r") as f:
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

pts = np.array(pts)
forward_reactions_serial = np.array(forward_reactions_serial)
forward_reactions_no_chunk = np.array(forward_reactions_no_chunk) 
forward_reactions_20 = np.array(forward_reactions_20) 
forward_reactions_100 = np.array(forward_reactions_100)

plt.plot(pts, forward_reactions_serial, 'k')
plt.show()
