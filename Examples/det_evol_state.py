import matplotlib.pyplot as plt
import numpy as np

filein = open("det_evol.txt")
tmin = float(filein.readline())
tmax = float(filein.readline())
dt = float(filein.readline())
npoints = int((tmax-tmin)/dt)

t = np.arange(tmin, tmax, dt)
psi = np.zeros(npoints)
phi = np.zeros(npoints)
nphi = np.zeros(npoints)
Phi = np.zeros(npoints)

for i in range(npoints):
  psi[i] = float(filein.readline())
  phi[i] = float(filein.readline())
  nphi[i] = float(filein.readline())
  Phi[i] = float(filein.readline())

fig, ax = plt.subplots(1,1)

ax.plot(t,psi, color = "red", label = r'$\psi$')
ax.plot(t,phi, color = "blue", label = r'$\phi$')
#ax.plot(t,Phi, color = "green", label = r'$\Phi$')
ax.set_ylabel(r'$tr[\rho\sigma_z]$')
ax.set_xlabel(r'$t$')

axx = ax.twinx()
axx.plot(t, nphi, '--', color = "green", label = r'$||\Phi||^2$')
axx.set_ylabel(r'$||\Phi||^2$')

ax.legend(loc = "lower right")
axx.legend(loc = "lower left")



plt.show()