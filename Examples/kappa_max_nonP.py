import matplotlib.pyplot as plt
import numpy as np

fig,ax = plt.subplots()

filein = open("kappa_theta_npoints.txt")
npoints = int(filein.readline())

filein = open("kappa_theta.txt")
kappa = np.zeros(npoints)
theta = np.zeros(npoints)
for i in range(npoints):
    pair = filein.readline()
    kappa[i] = float(pair.split()[0])
    theta[i] = float(pair.split()[1])


ax.scatter(kappa,theta,s=.01)
ax.set_xlabel(r'$\kappa$')
ax.set_ylabel(r'$\theta$')
ax.set_title(r'$|\psi_0> = \cos\theta|g> + \sin\theta|e>$')

plt.show()