import matplotlib.pyplot as plt
import numpy as np

filein = open("params.txt")
Ncopies = int(filein.readline())
Nensemble = int(filein.readline())
ti = float(filein.readline())
tf = float (filein.readline())
dt = float(filein.readline())
dimH = int(filein.readline())
Npoints = int((tf-ti)/dt)

t = np.arange(ti,tf,dt)

#observable = np.zeros((Nensemble, Npoints))
#exact = np.zeros(Npoints)
avg_obs = np.zeros(Npoints)
#filein = open("observables.txt")
#exact = open("analytic.txt")
f_avg = open("average.txt")
for i in range(Npoints):
    #noJumps[i] = f_noJ.readline()
    #exact[i] = f_exact.readline()
    avg_obs[i] = f_avg.readline()
    #j = 0
    #line = filein.readline()
    #for x in line.split():
        #observables[j,i] = x
        #j+=1
#for i in range(4):
    #plt.plot(t, observables[i,:], alpha=.1)
#plt.plot(t,exact,color='black', label="Exact", linewidth=1)
plt.plot(t,avg_obs,'g-', label="Average")

plt.legend(loc="lower right")
plt.ylabel(r'$Re(\rho_{12})$')
plt.xlabel(r'$t$')
plt.title('Driven')

plt.savefig("Eternally_nm.png")

plt.show()