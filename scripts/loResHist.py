import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


fig = plt.figure()
ax = fig.add_subplot(111)



# read in the data from zssa22a.current
data = np.genfromtxt("../spec/mask_design/z.ssa22a.current", skip_header=1)

print(data)
# loop through all of the rows in the data and add the redshifts to a list
zs = np.array([])
for row in data:
    zs = np.append(zs, float(row[8]))

print(np.where(zs>0))
N = len(np.where(zs > 0)[0])

    
plt.text(0.1, 0.9,"$N={}$".format(N), horizontalalignment='center', verticalalignment='center',transform = ax.transAxes, fontsize=16)

plt.hist(zs, range=(2, 3.7), bins=40, color='black')
plt.xlabel("$z$", fontsize=16)
plt.ylabel("$N$", fontsize=16)
plt.savefig("loResHistogram_LBG.png", dpi=400)
