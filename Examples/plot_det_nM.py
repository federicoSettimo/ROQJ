import matplotlib.pyplot as plt
import numpy as np

filein = open("det_nM_np.txt")
npoints = int(filein.readline())

Ndet = np.zeros(npoints)
Ntot = np.zeros(npoints)
kappa = np.zeros(npoints)

filein = open("det_nM.txt")
for i in range(npoints):
  line = filein.readline()
  kappa[i] = float(line.split()[0])
  Ntot[i] = float(line.split()[1])
  Ndet[i] = float(line.split()[2])

plt.plot(kappa, Ntot, label = "N tot")
plt.plot(kappa, Ndet, label = "N det")
plt.xlabel(r'$\kappa$')
plt.ylabel(r'$N$')
plt.legend(loc = "upper left")

plt.show()