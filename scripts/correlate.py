#Code to calculate the correlation function of galaxies
import numpy as np
from cosmology import *
import matplotlib.pyplot as plt
import sys

#Constants and conversions
cmToMpc = 3.241e-25
cmToGpc = 3.241e-28



#define class to hold galaxy data
class galaxy():
    def __init__(self, name, RA, dec, z):
        self.name = name
        self.RA = RA
        self.dec = dec
        self.z = z



#this will read in the galaxy data and return an dictionary of galaxy
# class objects
def galaxy_read():
    galaxies={} 
    #read in positional information 
    posData = np.genfromtxt("./masterPositions.dat", dtype=str) 
 
    #read in redshift data 
    zData = np.genfromtxt("table.dat", dtype=str, skip_header=1) 
 
    #loop through each galaxy that we have redshift information for and create a new     galaxy object 
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
                galaxies[zName] = galaxy(zName, float(posData[posIndex,1][0]), float(    posData[posIndex,2][0]), z) 
    print("Done adding galaxy data to {} galaxies.".format(len(galaxies))) 
    return galaxies 



#this function will calculate the 3d correlation function
# input the galaxy dictionary from galRead function
def three_dim_correlate(galaxies):
    #first loop through all of the galaxies and create a 3d array of positions
    positions = np.array([])
    for galKey in sorted(galaxies):
        galaxy = galaxies[galKey]
        # Define some variables for the extent of the cluster
        # these are roughly the middle of the cluster
        midRA = 334.37
        midDec =0.2425
        # this is the middle of the cluster in z space
        peakLimit = 3.078
        # each row will have the format [x, y, z]
        #  where x is the RA direction, y is the declination
        #  and the z component is the l.o.s. direction.
        x = cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.)
        y = cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.)
        z = cmToMpc*DC(galaxy.z)-cmToMpc*DC(peakLimit) 
       
        newline = [x, y, z]
        positions = np.append(positions, newline) 
    
    # reshape the array to be a 3xN array
    positions = np.reshape(positions, (-1, 3))

    toolbar_width=40
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))

    # define an array to hold the data galaxy correlations
    DD = np.array([])
    # define the distance range we are looking for
    dr = 0.4
    maxR = 40
    nsteps = maxR/dr
    # loop through all distances
    for ii in range(int(nsteps)):
        R = ii*dr
          
        #loop through each of the data galaxies
        Ntot = 0
        for thisGal in range(len(positions[:,0])):
            # create an empty array for the distances to this galaxy
            distances = np.array([])
            # loop through the other galaxies
            for otherGal in range(len(positions[:,0])):
                # make sure we are not counting the same galaxy
                if thisGal != otherGal:
                    # calculate the distance from thisgal to otherGal
                    distance = np.sum((positions[thisGal]-positions[otherGal])**2)
                    distances = np.append(distances, distance)

            # find the number that are within R-dr/2 , R+dr/2
            N = len(np.where(np.logical_and(distances>= R-dr/2., distances <= R+dr/2.))[0])
            Ntot += N         
        
        DD = np.append(DD, Ntot/2.)
        if ii%(nsteps/toolbar_width) == 0:
            sys.stdout.write("-")
            sys.stdout.flush()

    sys.stdout.write("\n")
    print(DD)
    plt.plot(np.linspace(0, maxR, nsteps),DD)
    plt.show()







if __name__=="__main__":
    print("----MAIN----")
    three_dim_correlate(galaxy_read())
