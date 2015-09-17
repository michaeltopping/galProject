import numpy as np
import matplotlib.pyplot as plt
import dataIO
from cosmology import *

from matplotlib.colors import LinearSegmentedColormap

#Constants and conversions
cmToMpc = 3.241e-25
cmToGpc = 3.241e-28

# plot all of the galaies and their redshifts
#  the galaxy shapes are separated by which peak they reside in
def plot_all_galaxies(categories, peak):
    # get the dictionary of all galaxies
    galaxies = dataIO.galaxy_init()
    #plot all of the positions in 2d
    #define colormap
    peakmin = 3.06
    peakmax = 3.11
    peakLimit = 3.078
    midPt = (peakLimit - peakmin) / (peakmax - peakmin)

    cdict = {'red':   ((0.0,  0.2, 0.0),
                       (midPt,  1.0, 1.0),
                       (1.0,  1.0, 1.0)),

             'green': ((0.0,  0.2, 0.0),
                       (midPt,  1.0, 1.0),
                       (1.0,  0.0, 0.0)),

             'blue':  ((0.0,  1.0, 1.0),
                       (midPt,  1.0, 1.0),
                       (1.0,  0.0, 0.0))}
    colorMap = LinearSegmentedColormap('BlueRed1', cdict)
    midRA = 334.37
    midDec =0.2425
    #loop through all galaxies
    for gal in sorted(galaxies):
        galaxy = galaxies[gal]
        #check which peak the galaxy lies in
        if galaxy.z > peakLimit:
            marker = 'o'
        else:
            marker = 'v'
        #make sure the galaxies are in the required range
        if peak == 'b':
            peakRange = [3.05, peakLimit]
        elif peak == 'r':
            peakRange = [peakLimit, 3.12]
        else:
            peakRange = [3.05, 3.12]
        if galaxy.z > peakRange[0] and galaxy.z < peakRange[1]:
            #plot the galaxy
            plt.scatter(cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.),
             cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.),
             vmin=cmToMpc*DC(peakmin)-cmToMpc*DC(peakLimit),
             vmax=cmToMpc*DC(peakmax)-cmToMpc*DC(peakLimit),
             c=cmToMpc*DC(galaxy.z)-cmToMpc*DC(peakLimit),
             cmap=colorMap, s=100, marker = marker)



    #plotting parameters
    plt.colorbar().set_label(r"$\rm Relative \ Comoving \ Distance \ [Mpc]$", fontsize=16)
    plotRange = 16
    plt.ylim([-plotRange/2.,plotRange/2.])
    plt.xlim([-plotRange/2.,plotRange/2.])
    plt.xlabel(r"$\rm Comoving \ Mpc$", fontsize=16)
    plt.ylabel(r'$\rm Comoving \ Mpc$', fontsize=16)


    #display/save figure
    plt.show()
   




















if __name__ == "__main__":
    print("MAIN")
    # 'peak' separates the symbol by what peak the galaxies are in
    # ' 
    plot_all_galaxies("peak",'d')
