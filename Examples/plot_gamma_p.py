# Plots both the trajectories and the functions
import matplotlib.pyplot as plt
import numpy as np
import sys

# Reading parameters
filein = open("params.txt")
Ncopies = int(filein.readline())
Nensemble = int(filein.readline())
ti = float(filein.readline())
tf = float (filein.readline())
dt = float(filein.readline())
print_traj = bool(filein.readline())
Ntraj = int(filein.readline())
dimH = int(filein.readline())
Npoints = int((tf-ti)/dt)

t = np.arange(ti,tf,dt)

# Reading trajectories and functions
trajectories = np.zeros((Ntraj, Npoints))
exact = np.zeros(Npoints)
avg_obs = np.zeros(Npoints)
err_obs = np.zeros(Npoints)
gp = np.zeros(Npoints)
gm = np.zeros(Npoints)
gz = np.zeros(Npoints)
if print_traj == True:
    filein = open("trajectories.txt")
f_exact = open("analytic.txt")
f_avg = open("average.txt")
f_err = open("error.txt")
f_func = open("functions.txt")
for i in range(Npoints):
    exact[i] = f_exact.readline()
    avg_obs[i] = f_avg.readline()
    err_obs[i] = f_err.readline()
    funcs = f_func.readline()
    gp[i] = float(funcs.split()[0])
    gm[i] = float(funcs.split()[1])
    gz[i] = float(funcs.split()[2])
    if print_traj == True:
        j = 0
        line = filein.readline()
        for x in line.split():
            trajectories[j,i] = x
            j+=1

# Plotting
fig, ax = plt.subplots(1,1, figsize = (13,4))

if print_traj == True:
    for i in range(Ntraj):
        ax.plot(t, trajectories[i,:], alpha=.1)
ax.plot(t,exact,color='black', label="Exact")
ax.errorbar(t,avg_obs,err_obs, marker='o', markersize=3, color='red', label="Average", errorevery=30, markevery=30, linewidth=0, elinewidth=1)

ax.legend(loc="upper left")
ax.set_xlabel(r'$t$')

# Inset: dephasing functions
ax2 = fig.add_axes([.71, .2, .18, .35])
ax2.plot(t, gp, label=r'$\gamma_+$')
ax2.plot(t, gm, label=r'$\gamma_-$')
ax2.plot(t, gz, label=r'$\gamma_z$')
ax2.legend(loc = "upper right")
ax2.set_xlabel(r'$t$')

if sys.argv.__len__() > 1:
    plt.suptitle(sys.argv[1])

if sys.argv.__len__() > 2:
    ax.set_ylabel(sys.argv[2])

if sys.argv.__len__() > 3:
    plt.savefig("Examples/"+sys.argv[3])

plt.show()