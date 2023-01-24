import matplotlib.pyplot as plt
import numpy as np

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

trajectories = np.zeros((Ntraj, Npoints))
exact = np.zeros(Npoints)
avg_obs = np.zeros(Npoints)
err_obs = np.zeros(Npoints)
filein = open("trajectories.txt")
f_exact = open("analytic.txt")
f_avg = open("average.txt")
f_err = open("error.txt")
for i in range(Npoints):
    exact[i] = f_exact.readline()
    avg_obs[i] = f_avg.readline()
    err_obs[i] = f_err.readline()
    j = 0
    line = filein.readline()
    for x in line.split():
        trajectories[j,i] = x
        j+=1
for i in range(Ntraj):
    plt.plot(t, trajectories[i,:], alpha=.1)
plt.plot(t,exact,color='black', label="Exact", linewidth=1)
plt.errorbar(t,avg_obs,err_obs,'g-', label="Average")
#plt.plot(t,avg_obs,'g-', label="Average")

plt.legend(loc="lower right")
plt.ylabel(r'$Re(\rho_{12})$')
plt.xlabel(r'$t$')
plt.title('Driven')

plt.savefig("Examples/Eternally_nm.png")

plt.show()