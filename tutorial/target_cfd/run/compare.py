import numpy as np
import matplotlib.pyplot as plt

dat1 = np.loadtxt("reactor_results.csv",skiprows=1,delimiter=',')
dat2 = np.loadtxt("chemFoam1.out",skiprows=1)
dat3 = np.loadtxt("chemFoam.out",skiprows=1)
plt.plot(dat1[:,0],dat1[:,1],'-k')
plt.plot(dat2[:,0],dat2[:,1],'--r')
plt.plot(dat3[:,0],dat3[:,1],'--g')
plt.show()
