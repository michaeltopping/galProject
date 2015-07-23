#ZPLOC
#Z Position LyA Correlator
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from matplotlib.colors import LinearSegmentedColormap


galaxies = {}
#define class to hold galaxy data
class galaxy():
	def __init__(self, name, RA, dec, z):
		self.name = name
		self.RA = RA
		self.dec = dec
		self.z = z


#read in positional information
posData = np.genfromtxt("./masterPositions.dat", dtype=str)

#read in redshift data
zData = np.genfromtxt("table.dat", dtype=str, skip_header=1)

#loop through each galaxy that we have redshift information for and create a new galaxy object
for ii in range(np.shape(zData)[0]):
	#the name of the galaxy
	zName = zData[ii][0] 
	#this is the systemic redshift
	z = float(zData[ii][1])
	#if there is no redshift information, z will be -9.999
	if z > 0:
		#the index in the position data array that contains this galaxy
		posIndex = np.where(zName == posData[:,0])[0]
		#check if the galaxy exists in the position array
		if posIndex:
			#create a new galaxy object in the galaxies dictionary
			galaxies[zName] = galaxy(zName, float(posData[posIndex,1][0]), float(posData[posIndex,2][0]), z)
print("Done adding galaxy data to {} galaxies.".format(len(galaxies)))




#-----Analysis-----

#-----Visualization-----
#plot all of the positions in 2d
#define colormap
cm = plt.get_cmap('seismic')
peakmin = 3.06
peakmax = 3.12
peakLimit = 3.08
midPt = (peakLimit - peakmin) / (peakmax - peakmin)

cdict = {'red':   ((0.0,  0.2, 0.0),
                   (midPt,  1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'green': ((0.0,  0.2, 0.0),
         		   (midPt,  1.0, 1.0),
                   (1.0,  0.0, 0.0)),

         'blue':  ((0.0,  1.0, 1.0),
                   (midPt,  1.0, 1.0),
                   (1.0,  0.2, 0.0))}
colorMap = LinearSegmentedColormap('BlueRed1', cdict)
#loop through all galaxies
for gal in galaxies:
	galaxy = galaxies[gal]
	#check which peak the galaxy lies in
	if galaxy.z > peakLimit:
		marker = 'o'
	else:
		marker = 'v'
	#make sure the galaxies are in the required range
	if galaxy.z > 3.05 and galaxy.z < 3.12:
		#plot the galaxy
		plt.scatter(galaxy.RA, galaxy.dec, vmin=peakmin, vmax=peakmax, c=galaxy.z, cmap=colorMap, s=100, marker = marker)

#plotting parameters
plt.colorbar().set_label("$z$", fontsize=16)
plt.xlim([334.315, 334.42])
plt.ylim([0.177,0.316])
plt.xlabel(r"$\rm Right \ Ascension \ [deg]$", fontsize=16)
plt.ylabel(r'$\rm Declination \ [deg]$', fontsize=16)
plt.savefig("scatter.png", dpi=400)
plt.show()
























