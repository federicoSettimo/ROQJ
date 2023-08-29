import matplotlib.pyplot as plt
import numpy as np

filein = open("target_state_unravelling.txt")
tmax = float(filein.readline())
dt = float(filein.readline())
npoints = int(tmax/dt)

t = np.arange(0,tmax,dt)
psix = np.zeros(npoints)
psiz = np.zeros(npoints)
phix = np.zeros(npoints)
phiz = np.zeros(npoints)
nphi = np.zeros(npoints)

for i in range(npoints):
  line = filein.readline()
  l = line.split()
  psix[i] = l[0]
  phix[i] = l[1]
  psiz[i] = l[2]
  phiz[i] = l[3]
  nphi[i] = l[4]

fig, ax = plt.subplots(1,1, figsize = (8,5))

ax.plot(t,psix, '-', color="red", label=r'$x$')
ax.plot(t,phix, '--', color="red")
ax.plot(t,psiz, '-', color="green", label=r'$z$')
ax.plot(t,phiz, '--', color="green")
ax.plot(t,nphi, '--', label = r'$||\Phi||^2$')
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$\psi$ solid, $\phi$ dashed')
ax.legend(loc = "lower right")
ax.set_title(r'$\gamma_{\pm}<0$')

plt.show()