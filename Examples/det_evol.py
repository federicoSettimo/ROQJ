import matplotlib.pyplot as plt
import numpy as np

# Deterministic evolution
filein = open("det_evol.txt")
nsteps = int(filein.readline())
nstates = int(filein.readline())
alpha0 = float(filein.readline())

t = np.zeros((nsteps, nstates))
alpha = np.zeros((nsteps, nstates))

for psi in range(nstates):
  for tt in range(nsteps):
    line = filein.readline()
    t[tt,psi] = float(line.split()[0])
    alpha[tt,psi] = float(line.split()[1])

# Alpha not ok
filein = open("alpha_no.txt")
alpha_no_min = float(filein.readline())
alpha_no_max = float(filein.readline())
tmin = float(filein.readline())
tmax = float(filein.readline())

# After-jump evolution
filein = open("det_evol_jump.txt")
tj = np.zeros((nsteps,nsteps))
psij = np.zeros((nsteps,nsteps))
for ttj in range(nsteps):
  for tt in range(nsteps):
    line = filein.readline()
    tj[tt,ttj] = float(line.split()[0])
    psij[tt,ttj] = float(line.split()[1])


# Plotting
fig, ax = plt.subplots(1,2, figsize = (11,3), sharey = True)

for psi in range(nstates):
  ax[0].plot(t[:,psi],alpha[:,psi], alpha=.2)

ax[0].axhline(alpha0,t[0,0],t[-1,0], color = "green", label = "Post jump state")
ax[0].set_ylabel(r'$\alpha = |<0|\psi>|$')
ax[0].set_xlabel(r'$t$')
ax[0].set_title("Deterministic evolution")
ax[0].axhspan(alpha_no_min, alpha_no_max, xmin = tmin/t[-1,0], xmax = tmax/t[-1,0], color = "red", alpha = .3, label = r'$\frac{1-\alpha^2}{\alpha^2} = \sqrt{\frac{\gamma_-}{\gamma_+}}$')
#ax[0].legend(loc = "lower right")
ax[0].text(-1,.97,r'$|0>$')
ax[0].text(-1,-.03,r'$|1>$')

for ttj in range(nsteps):
  ax[1].plot(tj[:,ttj], psij[:,ttj], alpha=.2)

ax[1].axhline(alpha0,t[0,0],t[-1,0], color = "green", label = "Post jump state")
ax[1].set_xlabel(r'$t$')
ax[1].set_title("Post-jump evolution")
ax[1].axhspan(alpha_no_min, alpha_no_max, xmin = tmin/t[-1,0], xmax = tmax/t[-1,0], color = "red", alpha = .3, label = r'$\frac{1-\alpha^2}{\alpha^2} = \sqrt{\frac{\gamma_-}{\gamma_+}}$')
ax[1].legend(loc = "lower left")

plt.show()