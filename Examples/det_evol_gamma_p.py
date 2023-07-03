import matplotlib.pyplot as plt
import numpy as np

# Deterministic evolution
filein = open("det_evol.txt")
feig = open("eigs.txt")
nsteps = int(filein.readline())
nstates = int(filein.readline())

t = np.zeros((nsteps, nstates))
z = np.zeros((nsteps, nstates))
x = np.zeros((nsteps, nstates))
eig = np.zeros((nsteps, nstates))

for psi in range(nstates):
  for tt in range(nsteps):
    line = filein.readline()
    t[tt,psi] = float(line.split()[0])
    z[tt,psi] = float(line.split()[1])
    x[tt,psi] = float(line.split()[2])
    eig[tt,psi] = float(feig.readline())

# After-jump evolution
filein = open("det_evol_jump.txt")
feig = open("eigs_jump.txt")

tj = np.zeros((nsteps,2*nsteps))
zj = np.zeros((nsteps,2*nsteps))
xj = np.zeros((nsteps,2*nsteps))
eigj = np.zeros((nsteps,2*nsteps))
for ttj in range(2*nsteps):
  for tt in range(nsteps):
    line = filein.readline()
    tj[tt,ttj] = float(line.split()[0])
    zj[tt,ttj] = float(line.split()[1])
    xj[tt,ttj] = float(line.split()[2])
    eigj[tt,ttj] = float(feig.readline())


# Plotting
fig, ax = plt.subplots(3,2, figsize = (11,6), sharey = False, sharex = True)

# Det evol
for psi in range(nstates):
  ax[0,0].plot(t[:,psi],z[:,psi], alpha=.2)
  ax[1,0].plot(t[:,psi],x[:,psi], alpha=.2)

ax[0,0].set_title("Deterministic evolution")
ax[0,0].set_ylabel(r'$tr[\rho\sigma_z]$')
ax[1,0].set_ylabel(r'$tr[\rho\sigma_x]$')

# Post-jump evol
for ttj in range(2*nsteps):
  if ttj < nsteps:
    ax[0,1].plot(tj[:,ttj], zj[:,ttj], alpha=.1, color = "green")
    ax[1,1].plot(tj[:,ttj], xj[:,ttj], alpha=.1, color = "green")
  else:
    ax[0,1].plot(tj[:,ttj], zj[:,ttj], alpha=.1, color = "blue")
    ax[1,1].plot(tj[:,ttj], xj[:,ttj], alpha=.1, color = "blue")

ax[0,1].set_title("Post-jump evolution")
ax[0,1].yaxis.set_ticklabels([])
ax[1,1].axhline(1,t[0,0],t[-1,0], color = "green", label = r'$|+>$')
ax[1,1].axhline(-1,t[0,0],t[-1,0], color = "blue", label = r'$|->$')
ax[1,1].legend(loc = "lower right")
ax[1,1].yaxis.set_ticklabels([])

# Eigenvalues
for psi in range(nstates):
  ax[2,0].plot(t[:,psi],eig[:,psi], alpha=.2)

ax[2,0].axhline(0,t[0,0],t[-1,0], color = "red")
ax[2,0].set_ylabel("Min eigenvalue")
ax[2,0].set_xlabel(r'$t$')
ax[2,0].set_ylim([-0.05,0.2])

# Eigenvalues post-jump
for ttj in range(2*nsteps):
  if ttj < nsteps:
    ax[2,1].plot(tj[:,ttj], eigj[:,ttj], alpha=.1, color = "green")
  else:
    ax[2,1].plot(tj[:,ttj], eigj[:,ttj], alpha=.1, color = "blue")

ax[2,1].axhline(0,t[0,0],t[-1,0], color = "red")
ax[2,1].set_xlabel(r'$t$')
ax[2,1].set_ylim([-0.05,0.2])
ax[2,1].yaxis.set_ticklabels([])

plt.tight_layout()
plt.show()