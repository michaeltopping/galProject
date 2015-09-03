import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from matplotlib.colors import LinearSegmentedColormap
from cosmology import *
from scipy import interpolate
import DBSCAN
from math import *


#Constants and conversions
cmToMpc = 3.241e-25
cmToGpc = 3.241e-28


#class to hold galaxy data
class galaxy():
    def __init__(self, name, RA, dec, z):
        self.name = name
        self.RA = RA
        self.dec = dec
        self.z = z
 
def dataInit():
    galaxies={}
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
                galaxies[zName] = galaxy(zName, float(posData[posIndex,1][0]), float(posData[posIndex,2]    [0]), z)
    print("Done adding galaxy data to {} galaxies.".format(len(galaxies)))
    return galaxies



def nearest(galaxies):

    midRA = 334.37
    midDec = 0.2425
    peakLimit = 3.077

    #put the positions of all galaxies into one array
    positions = np.array([])
    for gal in galaxies:
        galaxy = galaxies[gal]
        x = cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.)
        y = cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.)
        z = (cmToMpc*DC(galaxy.z)-cmToMpc*(DC(peakLimit)))
        positions = np.append(positions, [x, y, z])
    positions = np.reshape(positions, (-1, 3))
    #create the map to save the data
    nx = 100
    ny = 100

    xs = np.linspace(-8, 8, nx)
    yx = np.linspace(-8, 8, ny)
    map = np.zeros((nx, ny))

    #loop through each of the positions in the map
    for ix in range(nx):
        for iy in range(ny):
            xpositions = positions[:,0]
            ypositions = positions[:,1]
            zpositions = positions[:,2]
            distances = np.sqrt((xpositions-xs[ix])**2+(ypositions-yx[iy])**2+(zpositions-peakLimit)**2)
            print(min(distances)) 
           
            map[ix,iy] = 1/(min(distances)+10)
    plt.pcolor(map)
    plt.show()

if __name__=="__main__":
    galaxies = dataInit()
    nearest(galaxies)


