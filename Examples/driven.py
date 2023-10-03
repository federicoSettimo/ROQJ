import matplotlib.pyplot as plt
import numpy as np

filein = open("driven.txt")
Npoints = int(filein.readline())
tf = float(filein.readline())
dt = float(filein.readline())

t = np.arange(0,tf,dt)
psi = np.zeros(Npoints)
rho = np.zeros(Npoints)
exact = np.zeros(Npoints)
Npsi = np.zeros(Npoints)
Np = np.zeros(Npoints)
Nm = np.zeros(Npoints)
Ng = np.zeros(Npoints)
Ne = np.zeros(Npoints)
gp = np.zeros(Npoints)
gm = np.zeros(Npoints)
b = np.zeros(Npoints)

for i in range(Npoints):
  line = filein.readline()
  rho[i] = float(line.split()[0])
  exact[i] = float(line.split()[1])
  psi[i] = float(line.split()[2])

  line = filein.readline()
  Npsi[i] = float(line.split()[0])
  Ng[i] = float(line.split()[1])
  Ne[i] = float(line.split()[2])
  Np[i] = float(line.split()[3])
  Nm[i] = float(line.split()[4])

  line = filein.readline()
  gp[i] = float(line.split()[0])
  gm[i] = float(line.split()[1])
  b[i] = float(line.split()[2])

fig, ax = plt.subplots(1,3, figsize=(17,4))

ax[0].plot(t, rho, color = "red", label = r'$\rho$')
ax[0].plot(t, exact, color = "black")
ax[0].plot(t ,psi, '--', label = r'$\psi$', alpha = .5)
ax[0].set_ylabel(r'$tr[\rho\sigma_z]$')
ax[0].set_xlabel(r'$t$')
ax[0].legend(loc = "lower right")

ax[1].plot(t, Npsi, label = r'$\psi$')
ax[1].plot(t, Ng, label = r'$0$')
ax[1].plot(t, Ne, label = r'$1$')
ax[1].plot(t, Np, label = r'$+$')
ax[1].plot(t, Nm, label = r'$-$')
ax[1].set_ylabel(r'Populations')
ax[1].set_xlabel(r'$t$')
ax[1].legend(loc = "lower right")

ax[2].plot(t,gp, label = r'$\gamma_+$')
ax[2].plot(t,gm, label = r'$\gamma_-$')
ax[2].plot(t,b, label = r'$\beta$')
ax[2].set_xlabel(r'$t$')
ax[2].legend(loc = "lower right")

plt.show()