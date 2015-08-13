#ZPLOC
#Z Position LyA Correlator
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from matplotlib.colors import LinearSegmentedColormap
from cosmology import *
from scipy import interpolate
import DBSCAN


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
                galaxies[zName] = galaxy(zName, float(posData[posIndex,1][0]), float(posData[posIndex,2][0]), z)
    print("Done adding galaxy data to {} galaxies.".format(len(galaxies)))
    return galaxies



#-----Analysis-----

#-----Visualization-----i
def plot2d(galaxies,peak):
    #plot all of the positions in 2d
    #define colormap
    peakmin = 3.06
    peakmax = 3.11
    peakLimit = 3.077
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
    for gal in galaxies:
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
    plt.savefig("scatter.png", dpi=400)
    plt.show()



def plot3d(galaxies, movieBool):
    #Basic plotting vaiables
    peakmin = 3.06
    peakmax = 3.11
    peakLimit = 3.077
    midRA = 334.37
    midDec =0.2425
    
    #Create figure for 3d visualization
    fig3d = plt.figure()
    ax3d = fig3d.add_subplot(111, projection='3d')

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
            #choose the color of the point based on which peak its in
            if galaxy.z < peakLimit:
                color = 'blue'
                marker = 'v'
            if galaxy.z > peakLimit:
                color = 'red'
                marker = 'o'
            #plot the galaxy
            ax3d.scatter(cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.),
                 cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.),
                 (cmToMpc*DC(galaxy.z)-cmToMpc*(DC(peakLimit))), color=color, 
                 marker=marker, edgecolors='black')
    #plotting parameters
    ax3d.set_zlim3d(-30, 30)
    ax3d.set_xlim3d(-8, 8)
    ax3d.set_ylim3d(-8, 8)
    ax3d.set_xlabel(r"$\rm Comoving \ Mpc$")
    ax3d.set_ylabel(r"$\rm Comoving \ Mpc$")
    ax3d.set_zlabel(r"$\rm Comoving \ Mpc$")
    if movieBool:
        #Animate the azimuth for the 3d plot
        for ii in range(0,360,1):
            ax3d.view_init(elev=10., azim=ii)
            plt.savefig("./movie/{}.png".format(ii), dpi=400)
            print("Finished saving frame {}.".format(ii))
    else:   
        plt.show()





#calculate the column density and plot it
def colDensity(galaxies, peak):
    #basic range parameters
    peakmin = 3.06
    peakmax = 3.11
    peakLimit = 3.077
    midRA = 334.37
    midDec =0.2425

    #dependeing on input, set the limits to enclose one or both of the peaks
    if peak == 'b':
        peakRange = [3.05, peakLimit]
    elif peak == 'r':
        peakRange = [peakLimit, 3.12]
    else:
        peakRange = [3.05, 3.12]
    #blank positions array 
    positions = np.array([])
    #loop through the dictionary of galaxies
    for gal in galaxies:
        galaxy = galaxies[gal] 

        #only look for galaxies that are in the requested peak
        if galaxy.z > peakRange[0] and galaxy.z < peakRange[1]:

            #create an array of x and y values to the galaxies 
            x = cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.)
            y = cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.)
            positions = np.append(positions, np.array([x, y]) )

    positions = positions.reshape(-1,2)
    #parameters for the size of the heatmap
    N =200 
    minx = -8.
    miny = -8.
    maxx = 8.
    maxy = 8.
    dx = (maxx-minx)/N
    dy = (maxy-miny)/N
    Rs = 20.
    #create the blank heat map
    density = np.zeros(shape=(N, N))
    
    #Loop through the x positions
    for ix in range(N):
        #loop through the y positions
        for iy in range(N):
            distances = np.array([])
            #calculate a list of distances
            distances = np.sqrt((positions[:,0]-(minx+ix*dx))**2+(positions[:,1]-(miny+iy*dy))**2)
            #Count the number of galaxies that are within 2 Comoving Mpc
#            density[ix][iy] = len(np.where(distances<2)[0])
            #divide the number of galaxies in the peak by the sum of the distances
            density[ix][iy] = np.sum(np.shape(distances)[0]/(distances/Rs*(1+distances/Rs)**2))

   #Create the plot
    std = np.std(density)
    plt.pcolor(density, vmin=0, vmax=5*std)
    plt.xlabel(r"$\rm Comoving \ Mpc$", fontsize=16)
    plt.ylabel(r"$\rm Comoving \ Mpc$", fontsize=16)
    plt.colorbar(label="{}".format(peak))
    plt.show()





#find the number of clusters in the sample
def clustering(galaxies, eps, minPts):
    midRA = 334.37
    midDec =0.2425
    peakLimit = 3.077
    #create the list of points
    D = np.array([])
    for gal in galaxies:
        galaxy = galaxies[gal]
        x = cmToMpc*DC(galaxy.z)*(galaxy.RA-midRA)*(np.pi/180.)
        y = cmToMpc*DC(galaxy.z)*(galaxy.dec-midDec)*(np.pi/180.)
        z = (cmToMpc*DC(galaxy.z)-cmToMpc*(DC(peakLimit)))
        D = np.append(D,[x, y, z])
    #reshape to create a 2xN array
    D = D.reshape((-1,3))

    clusters, noise = DBSCAN.DBSCAN(D, eps, minPts)
    return clusters, noise


#plot clusters and noise particles
def plotClusters(clusters, noise):
    

    colors = ['k', 'b', 'r', 'y', 'c']
     #loop through and plot the points
    for ii in range(len(clusters)):
        for P in clusters[ii].points:
            plt.plot(P.x, P.y, '{}^'.format(colors[ii%4]))
 
    #plot the noise points
    for P in noise:
        plt.plot(P.x, P.y, 'k.')
 
    plt.ylim([-8, 8])
    plt.xlim([-8, 8])
    plt.show()




def clusterParamSearch(galaxies):
    eArr = np.linspace(0, 20, 1000)
    nclusters = np.array([])    
    nnoise = np.array([])

    #loop through each of the parameters
    for e in eArr:
        cs, ns = clustering(galaxies, e, 3)
        nclusters = np.append(nclusters, len(cs))
        nnoise = np.append(nnoise, len(ns))

    plt.plot(eArr, nclusters, 'k', linewidth=2)
    plt.plot(eArr, nnoise, 'b', linewidth=2)
    plt.ylabel("Number of Clusters")
    plt.xlabel("eps")
    plt.show()



if __name__=="__main__":

    galaxies = dataInit()    
    while True:
        print("----Commands----")
        print("q: quit")
        print("2: plot 2d figure")
        print("3: plot 3d figure")
        print("4: calculate clustering")
        print("5: search clustering parameter space")
        char = input("Make Selection:")
        if char=='q':
            break
        if char=='2':
            print("Scatter plot or heat map? (s/h)")
            plottype=input(":")
            print("Which peak to plot?")
            print("  r: red peak")
            print("  b: blue peak")
            print("  d: both peaks")
            peak = input(":")
            if plottype == "s":
                plot2d(galaxies, peak)
            if plottype == "h":
                colDensity(galaxies, peak)
        if char=='3':
            movieBool = input("Create Movie? (y/n)") is 'y'
            plot3d(galaxies, movieBool)
        if char=='4':
            clustering(galaxies, 1, 3)
        if char=='5':
            clusterParamSearch(galaxies)







