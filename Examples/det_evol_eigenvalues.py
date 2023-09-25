import matplotlib.pyplot as plt
import numpy as np

filein = open("det_evol.txt")
tmin = float(filein.readline())
tmax = float(filein.readline())
dt = float(filein.readline())
npoints = int((tmax-tmin)/dt)

t = np.arange(tmin, tmax, dt)
psi = np.zeros(npoints)
eig1 = np.zeros(npoints)
eig2 = np.zeros(npoints)

for i in range(npoints):
  psi[i] = float(filein.readline())
  eig1[i] = float(filein.readline())
  eig2[i] = float(filein.readline())

fig, ax = plt.subplots(1,1)

ax.plot(t,eig1)
ax.plot(t,eig2)
ax.plot(t,np.zeros(npoints), color = "red", alpha = .5)
ax.set_ylabel(r'Eigenvalues (solid)')
ax.set_xlabel(r'$t$')

axx = ax.twinx()
axx.plot(t, psi, '--', color = "green", alpha = .3)
axx.set_ylabel(r'$tr[\rho\sigma_z]$ (dashed)')


plt.show()