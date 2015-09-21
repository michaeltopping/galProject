import numpy as np
import matplotlib.pyplot as plt
import dataIO
from cosmology import *
from collections import OrderedDict
from matplotlib import rc

from matplotlib.colors import LinearSegmentedColormap

# change typesetting for plots
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


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
        # get the marker for each galaxy
        marker, label = galaxy_marker(galaxy, categories)
        #make sure the galaxies are in the required range
        if peak == 'b':
            peakRange = [3.05, peakLimit]
        elif peak == 'r':
            peakRange = [peakLimit, 3.12]
        else:
            peakRange = [3.05, 3.12]
        if galaxy.z > peakRange[0] and galaxy.z < peakRange[1]:
            #plot the galaxy
            # make sure that the galaxy is not a quasar
            if 'q' not in galaxy.LyAtype:
                plt.scatter(cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.),
                 cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.),
                 vmin=cmToMpc*DC(peakmin)-cmToMpc*DC(peakLimit),
                 vmax=cmToMpc*DC(peakmax)-cmToMpc*DC(peakLimit),
                 c=cmToMpc*DC(galaxy.z)-cmToMpc*DC(peakLimit),
                 cmap=colorMap, s=100, marker = marker, label=label)

        # do this if the galaxy is a QSO
        if 'q' in galaxy.LyAtype:
            # galaxy is a quasar
             plt.scatter(cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.),
             cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.),
             vmin=cmToMpc*DC(peakmin)-cmToMpc*DC(peakLimit),
             vmax=cmToMpc*DC(peakmax)-cmToMpc*DC(peakLimit),
             c = 'black',
             cmap=colorMap, s=100, marker = '*', label=label)

           

    
    #plotting parameters
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper left')
    plt.colorbar().set_label(r"$\rm Relative \ Comoving \ Distance \ [Mpc]$", fontsize=16)
    plotRange = 18
    plt.ylim([-plotRange/2.,plotRange/2.])
    plt.xlim([-plotRange/2.,plotRange/2.])
    plt.xlabel(r"$\rm Comoving \ Mpc$", fontsize=16)
    plt.ylabel(r'$\rm Comoving \ Mpc$', fontsize=16)


    #display/save figure
    plt.savefig("scatter_"+categories+".png", dpi=400)
    plt.close()
   




# gets the galaxy marker
def galaxy_marker(galaxy, categories):
    marker = 'o'
    label="Missing"
    if categories == 'peak':
        peakLimit = 3.078
        # if the galaxies will be separated out by peak
        if galaxy.z < peakLimit:
            # the galaxy is in the blue peak
            marker = '^'
            label="Blue Peak"
        else:
            # the galaxy is in the red peak
            marker = 'o'
            label="Red Peak"
    elif categories == 'LAE_LBG':
        # if the galaxies will be separated by their LAE.LBG status
        if galaxy.name.startswith("C") or galaxy.name.startswith("D") or galaxy.name.startswith("M"):
            # the galaxy is an LBG
            marker = '^'
            label="LBG"
        elif galaxy.name.startswith("NB") or galaxy.name.startswith("lae") or galaxy.name[0].isdigit():
            # the galaxy is an LAE  
            marker = 'o'
            label="LAE"
    elif categories == 'LyA_peakType':
        if galaxy.LyAtype == 'd':
            marker = [[0,-1], [-.8,-1], [-.9, 1], [.9,1], [.8,-1],[0,-1],[0,1]]
            label="Double peaked LyA"
        else:
            marker = [[.5, -1], [0,-1], [-.5, -1], [-.4, 1], [.4, 1], [.5, -1]]
            label=r"Single peaked Ly$\alpha$"
    elif categories == 'lineType':
        if ( galaxy.emZ > 0 ) and ( galaxy.absZ > 0 ):
            # galaxy has both emission and absorption lines
            marker = [[-1.06, 0], [-0.35, 1], [0.35, 0], [1.06, 0], 
                        [0.35, -1], [-0.35, 0], [-1.06,0],[1.06, 0]]
            label=r"$z_{em} \ \rm and \ z_{abs}$"
        elif ( galaxy.emZ > 0 ):
            # galaxy has only emission lines
            marker = '^'
            label=r"$z_{em} \ \rm only$"
        elif ( galaxy.absZ > 0 ):
            # galaxy has only absorption lines
            marker = 'v'
            label=r"$z_{abs}$"
        else:
            print("Error, galaxy has no measured redshifts")
            

    else:
        print("Invalid category, exiting now.")
        exit() 


    return marker, label









# create a method that will plot the positions of galaxies from
#  yamada et al. 2012
def yamada_scatter():
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
    # file where the table is stored
    filename = "yamada_table.dat"
    # read in the data
    data = np.genfromtxt(filename, dtype=str)
    
    for row in data:
        RA = 15*( float(row[1][1:3])+float(row[1][3:5])/60.+\
                    float(row[1][5:])/3600. )
        dec = float(row[2][0:2])+float(row[2][2:4])/60.+\
                    float(row[2][4:])/3600.

        z = float(row[3])
        plt.scatter(cmToMpc*DC(z)*(RA-midRA)*(np.pi/180.),
         cmToMpc*DC(z)*(dec-midDec)*(np.pi/180.),
             vmin=cmToMpc*DC(peakmin)-cmToMpc*DC(peakLimit),
             vmax=cmToMpc*DC(peakmax)-cmToMpc*DC(peakLimit),
             c=cmToMpc*DC(z)-cmToMpc*DC(peakLimit),
             cmap=colorMap, s=100, marker = 's', label="Yamada et al. 2012")




if __name__ == "__main__":
    print("MAIN")
    # 'peak' separates the symbol by what peak the galaxies are in
    # 'LAE_LBG' separates out based on the galaxy type 
    # 'LyA_peakType' separates out the single and the double peaked LyA galaxis
    # 'lineType' will separate out galaxies that have zem, zabs, or both
    categoryList = ['peak', 'LAE_LBG', 'LyA_peakType', 'lineType']
    for category in categoryList:
        plot_all_galaxies(category, 'd')
#    yamada_scatter()
#    plot_all_galaxies('peak', 'd')
#    plt.savefig("scatter_yamada_color_small.png", dpi=400)
#    plt.show()
