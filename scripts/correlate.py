#Code to calculate the correlation function of galaxies
import numpy as np
from cosmology import *
import matplotlib.pyplot as plt
import sys
import random

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
        peakMin = 3.06  
        peakMax = 3.11
        # each row will have the format [x, y, z]
        #  where x is the RA direction, y is the declination
        #  and the z component is the l.o.s. direction.
        x = cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.)
        y = cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.)
        z = cmToMpc*DC(galaxy.z)-cmToMpc*DC(peakLimit) 
       
        newline = [x, y, z]
        if galaxy.z < peakMax and galaxy.z > peakMin:
            positions = np.append(positions, newline) 
    
    # reshape the array to be a 3xN array
    positions = np.reshape(positions, (-1, 3))


    # define an array to hold the data galaxy correlations
    DD = np.array([])
    # define the distance range we are looking for
    dr = 1.0
    maxR = 60
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
        sys.stdout.write("\rDD:  {}%".format(100.*(ii+1)/nsteps))

    sys.stdout.write("\n")


    # now create a correlation function for a random sample of galaxies
    DR = np.array([])
    nGalaxies = 500
    xmin, xmax = min(positions[:,0]), max(positions[:,0])
    ymin, ymax = min(positions[:,1]), max(positions[:,1])
    zmin, zmax = min(positions[:,2]), max(positions[:,2]) 

    # create a positions list for the random sample of galaxies
    randPositions = np.array([])
    for ii in range(nGalaxies):
        x = random.uniform(xmin, xmax)
        y = random.uniform(ymin, ymax)
        z = random.uniform(zmin, zmax)
        newline = [x, y, z]
        randPositions = np.append(randPositions, newline)
    randPositions = np.reshape(randPositions, (-1, 3))


    # now do this for all of the random positions
    # loop through all distances
    for ii in range(int(nsteps)):
        R = ii*dr
          
        #loop through each of the data galaxies
        Ntot = 0
        for thisGal in range(len(randPositions[:,0])):
            # create an empty array for the distances to this galaxy
            distances = np.array([])
            # loop through the other galaxies
            for otherGal in range(len(positions[:,0])):
                # make sure we are not counting the same galaxy
                if thisGal != otherGal:
                    # calculate the distance from thisgal to otherGal
                    distance = np.sqrt(np.sum((randPositions[thisGal]-positions[otherGal])**2))
                    distances = np.append(distances, distance)

            # find the number that are within R-dr/2 , R+dr/2
            N = len(np.where(np.logical_and(distances>= R-dr/2., distances <= R+dr/2.))[0])
            Ntot += N         
        
        DR = np.append(DR, Ntot/2.)
        sys.stdout.write("\rDR:  {}%".format(100.*(ii+1)/nsteps))
    sys.stdout.write("\n")

    xi = nGalaxies/len(positions[:,0])*DD/DR-1
    return xi






# this function will calculate the two-point correlation function
#  but will differentiate between radial and transverse directions.
def redshift_correlation(galaxies):


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
        peakMin = 3.06  
        peakMax = 3.11
        # each row will have the format [x, y, z]
        #  where x is the RA direction, y is the declination
        #  and the z component is the l.o.s. direction.
        x = cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.)
        y = cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.)
        z = cmToMpc*DC(galaxy.z)-cmToMpc*DC(peakLimit) 
       
        newline = [x, y, z]
        if galaxy.z < peakMax and galaxy.z > peakMin:
            positions = np.append(positions, newline) 
    
    # reshape the array to be a 3xN array
    positions = np.reshape(positions, (-1, 3))

    # now create a correlation function for a random sample of galaxies
    DR = np.array([])
    nGalaxies = 500
    xmin, xmax = min(positions[:,0]), max(positions[:,0])
    ymin, ymax = min(positions[:,1]), max(positions[:,1])
    zmin, zmax = min(positions[:,2]), max(positions[:,2]) 

    # create a positions list for the random sample of galaxies
    randPositions = np.array([])
    for ii in range(nGalaxies):
        x = random.uniform(xmin, xmax)
        y = random.uniform(ymin, ymax)
        z = random.uniform(zmin, zmax)
        newline = [x, y, z]
        randPositions = np.append(randPositions, newline)
    randPositions = np.reshape(randPositions, (-1, 3))



    # set the boundary limits
    rTransMax = 40
    rRadMax = 40
    dr = 0.5
    nsteps = rTransMax*rRadMax/dr
    ii = 0 

    # create variables to hold results
    DD = np.array([])
    DR = np.array([])


    # first loop through the transverse direction
    for rTrans in np.arange(0, rTransMax, dr):

        # now loop through each of the radial direction bins
        for rRad in np.arange(0, rRadMax, dr):
            ii += 1
            #loop through each of the data galaxies
            NRadtot = 0
            NTranstot = 0
            Ntot = 0
            for thisGal in range(len(positions[:,0])):
                # create an empty array for the distances to this galaxy
                radDistances = np.array([])
                transDistances = np.array([])
                # loop through the other galaxies
                for otherGal in range(len(positions[:,0])):
                    # make sure we are not counting the same galaxy
                    if thisGal != otherGal:
                        # calculate the distance from thisgal to otherGal
                        radDistance = np.sqrt((positions[thisGal][2]-positions[otherGal][0])**2)
                        transDistance = np.sqrt((positions[thisGal][0] -
                                         positions[otherGal][0])**2 + 
                                         (positions[thisGal][1] - 
                                         positions[otherGal][1])**2)
                        radDistances = np.append(radDistances, radDistance)
                        transDistances = np.append(transDistances, transDistance)

                # find the number that are within R-dr/2 , R+dr/2
                N = len(np.where(np.logical_and(np.logical_and( 
                    radDistances>= rRad-dr/2., radDistances <= rRad+dr/2. ), 
                    np.logical_and( transDistances>= rTrans-dr/2.,
                    transDistances <= rTrans+dr/2. ))[0]))
                Ntot += N         
            # store the result
            DD = np.append(DD, Ntot/2.)
 
        

            # now compute the correlation between data and random galaxies
            #loop through each of the data galaxies
            NRadtot = 0
            NTranstot = 0
            Ntot = 0
            for thisGal in range(len(randPositions[:,0])):
                # create an empty array for the distances to this galaxy
                radDistances = np.array([])
                transDistances = np.array([])
                # loop through the other galaxies
                for otherGal in range(len(positions[:,0])):
                    # make sure we are not counting the same galaxy
                    if thisGal != otherGal:
                        # calculate the distance from thisgal to otherGal
                        radDistance = np.sqrt((randPositions[thisGal][2]-positions[otherGal][0])**2)
                        transDistance = np.sqrt((randPositions[thisGal][0] -
                                         positions[otherGal][0])**2 + 
                                         (randPositions[thisGal][1] - 
                                         positions[otherGal][1])**2)
                        radDistances = np.append(radDistances, radDistance)
                        transDistances = np.append(transDistances, transDistance)

                # find the number that are within R-dr/2 , R+dr/2
                N = len(np.where(np.logical_and(np.logical_and( 
                    radDistances>= rRad-dr/2., radDistances <= rRad+dr/2. ), 
                    np.logical_and( transDistances>= rTrans-dr/2.,
                    transDistances <= rTrans+dr/2. ))[0]))
                Ntot += N         
                # store the result
                DR = np.append(DR, Ntot/2.)

            sys.stdout.write("\rDD:  {:.2f}%".format(100.*(ii+1)/nsteps))
    DD = np.reshape(DD, (rTransMax/dr, rRadMax/dr))
    DR = np.reshape(DR, (rTransMax/dr, rRadMax/dr))

    xi = nGalaxies/len(positions[:,0])*DD/DR-1
    return xi


if __name__=="__main__":
    print("----MAIN----")
    galaxies = galaxy_read()
#    totalXi = three_dim_correlate(galaxies)
    totalXi = redshift_correlation(galaxies)
    print(totalXi)
#    nSamples = 10
#    for ii in range(nSamples-1):
#        print("Begin sample number: {}".format(ii+1))
#        totalXi += three_dim_correlate(galaxies)
#
#    print(totalXi/nSamples)
#    print(len(totalXi))
#    plt.plot(np.linspace(0, 40, len(totalXi)), totalXi/nSamples, 'ko')
#    plt.xlabel('Comoving Mpc')
#    plt.ylabel(r'$\xi$')
#    plt.show()


