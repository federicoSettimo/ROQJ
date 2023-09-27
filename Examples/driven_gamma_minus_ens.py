import matplotlib.pyplot as plt
import numpy as np

filein = open("driven_gamma_minus_ens.txt")
Npoints = int(filein.readline())
tf = float(filein.readline())
dt = float(filein.readline())

t = np.arange(0,tf,dt)
psi = np.zeros(Npoints)
phi = np.zeros(Npoints)
phi_perp = np.zeros(Npoints)
rho = np.zeros(Npoints)
exact = np.zeros(Npoints)

for i in range(Npoints):
  line = filein.readline()
  rho[i] = float(line.split()[0])
  psi[i] = float(line.split()[1])
  phi[i] = float(line.split()[2])
  phi_perp[i] = float(line.split()[3])
  exact[i] = float(line.split()[4])

plt.plot(t, rho, color = "red", label = r'$\rho$')
plt.plot(t, exact, color = "black")
plt.plot(t ,psi, '--', label = r'$\psi$', alpha = .5)
plt.plot(t ,phi, '--', label = r'$\varphi$', alpha = .5)
plt.plot(t ,phi_perp, '--', label = r'$\varphi_{\perp}$', alpha = .5)

plt.ylabel(r'$tr[\rho\sigma_z]$')
plt.xlabel(r'$t$')
plt.legend(loc = "lower right")

plt.show()