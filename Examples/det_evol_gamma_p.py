import matplotlib.pyplot as plt
import numpy as np

# Deterministic evolution
filein = open("det_evol.txt")
feig = open("eigs.txt")
nsteps = int(filein.readline())
nstates = int(filein.readline())

t = np.zeros((nsteps, nstates))
alpha = np.zeros((nsteps, nstates))
eig = np.zeros((nsteps, nstates))

for psi in range(nstates):
  for tt in range(nsteps):
    line = filein.readline()
    t[tt,psi] = float(line.split()[0])
    alpha[tt,psi] = float(line.split()[1])
    eig[tt,psi] = float(feig.readline())

# After-jump evolution
filein = open("det_evol_jump.txt")
feig = open("eigs_jump.txt")

tj = np.zeros((nsteps,2*nsteps))
psij = np.zeros((nsteps,2*nsteps))
eigj = np.zeros((nsteps,2*nsteps))
for ttj in range(2*nsteps):
  for tt in range(nsteps):
    line = filein.readline()
    tj[tt,ttj] = float(line.split()[0])
    psij[tt,ttj] = float(line.split()[1])
    eigj[tt,ttj] = float(feig.readline())


# Plotting
fig, ax = plt.subplots(2,2, figsize = (11,3), sharey = False, sharex = False)

# Det evol
for psi in range(nstates):
  ax[0,0].plot(t[:,psi],alpha[:,psi], alpha=.2)

ax[0,0].set_ylabel(r'$tr[\rho\sigma_z]$')
#ax[0,0].set_xlabel(r'$t$')
ax[0,0].set_title("Deterministic evolution")

# Post-jump evol
for ttj in range(2*nsteps):
  if ttj < nsteps:
    ax[0,1].plot(tj[:,ttj], psij[:,ttj], alpha=.1, color = "green")
  else:
    ax[0,1].plot(tj[:,ttj], psij[:,ttj], alpha=.1, color = "blue")

ax[0,1].axhline(0,t[0,0],t[-1,0], color = "green", label = r"$|+>$")
ax[0,1].axhline(0,t[0,0],t[-1,0], color = "blue", label = r"$|->$")
#ax[0,1].set_xlabel(r'$t$')
ax[0,1].set_title("Post-jump evolution")
ax[0,1].legend(loc = "lower right")

# Eigenvalues
for psi in range(nstates):
  ax[1,0].plot(t[:,psi],eig[:,psi], alpha=.2)

ax[1,0].axhline(0,t[0,0],t[-1,0], color = "red")
ax[1,0].set_ylabel("Min eigenvalue")
#ax[1,0].set_xlabel(r'$t$')

# Eigenvalues post-jump
for ttj in range(2*nsteps):
  if ttj < nsteps:
    ax[1,1].plot(tj[:,ttj], eigj[:,ttj], alpha=.1, color = "green")
  else:
    ax[1,1].plot(tj[:,ttj], eigj[:,ttj], alpha=.1, color = "blue")

ax[1,1].axhline(0,t[0,0],t[-1,0], color = "red")
ax[1,1].set_xlabel(r'$t$')

plt.show()