import numpy as np
import random
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from math import *
from cosmology import *
import scipy
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap 
from matplotlib import cm
import operator
from matplotlib import rc


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


#Constants and conversions
cmToMpc = 3.241e-25
cmToGpc = 3.241e-28




#define two gaussians
def twoGauss(x, a1, c1, s1, a2, c2, s2):
    if s1 ==0 or s2 == 0:
        return 99999.9
    else:
        return a1*np.exp(-(x-c1)**2/(2*s1**2))+a2*np.exp(-(x-c2)**2/(2*s2**2))

#define Gauss
def Gauss(x, a1, c1, s1):
    if s1 ==0:
        return 99999.9
    else:
        return a1*np.exp(-(x-c1)**2/(2*s1**2))




# split the results of a query into files for each resultant halo
def split_query(filename):
    # read in the entire sql query result
    data = np.genfromtxt(filename, dtype=str, skip_header=55, delimiter=',')
    
    # find the root of the filename
    fileroot, fileext = os.path.splitext(filename)
    nCluster = 0
    outFilename = fileroot+"_"+str(nCluster)+".dat"
    outFile = open(outFilename, 'w')
    
    # keep track of which current descendent halo we are looking at
    currentDesHalo = data[0][0]

    # loop through the rows  
    for row in data:
        if not row[0]==currentDesHalo:
            currentDesHalo = row[0]
            nCluster += 1
            outFilename = fileroot+"_"+str(nCluster)+".dat"
            outFile.close()
            outFile = open(outFilename, 'w')

        outFile.write("{}   {}    {}    {}    {}    {}    {}    {}    {}\n".format(row[0], row[19], row[20],row[21], row[23], row[24], row[25], row[15], random.random()))
        
    


# read in the halo data from one file
def read_halos(folder, filename):
    # first read in the data from one of the files
    data = np.genfromtxt(filename, dtype=str)
    mp = 8.6e8


    positions = np.array([])
    velocities = np.array([])
    masses = np.array([])
    randoms = np.array([])
    for row in data:
        # get the line in the positions array that will be added
        newline = [float(row[1]), float(row[2]), float(row[3])]
        newVelLine = [float(row[4]), float(row[5]), float(row[6])]
        positions = np.append(positions, newline)
        velocities = np.append(velocities, newVelLine)
        masses = np.append(masses, float(row[7])*mp)
        randoms = np.append(randoms, float(row[8]))
        

    # reshape the array to be a 3xN array
    positions = np.reshape(positions, (-1,3))
    velocities = np.reshape(velocities, (-1, 3))
    
    # get the average value for each of the 3 dimensions
    avgx = np.average(positions[:,0]) 
    avgy = np.average(positions[:,1])
    avgz = np.average(positions[:,2])
   
    # now compute the positions relative to the average position
    relPositions = np.array([])
    for row in positions:
        relPositions = np.append(relPositions, [row[0]-avgx, row[1]-avgy, 
                                                row[2]-avgz])

    relPositions = np.reshape(relPositions, (-1,3))


    nTheta = 24
    nPhi = 24
    thetas = np.linspace(0, 2*np.pi, nTheta, endpoint=False)
    phis = np.linspace(0, 2*np.pi, nPhi, endpoint=False)

    # create the fs array, which will be thetas x phis
    fs = np.array([])    
    amplitudes = np.array([])
    widths = np.array([])

 
    vertList = get_icovertices(0, thetas, phis) 
    N = len(vertList)
    ii = 0
    totStartTime = time.time()
    for vert in vertList:
        starttime = time.time()
        # note that the parameters for single and double peaked fits
        #  are switched here because I was to lazy to switch all following references
        #  to them. :(
        f, param2, param1 = compute_histogram(np.array(rotate_galaxies(relPositions, vert).T), 
                            np.array(rotate_galaxies(velocities, vert).T), masses, vert, randoms)

#        scatter(np.array(rotate_galaxies(relPositions, vert).T),
#                            np.array(rotate_galaxies(velocities, vert).T), masses, vert, randoms)

        fs = np.append(fs, f)
        amplitudes = np.append(amplitudes, min(abs(param2[3]/param2[0]), abs(param2[0]/param2[3])))
        widths = np.append(widths, min(abs(param2[2]/param2[5]), abs(param2[5]/param2[2])))
        print("Finished raytrace: {} in {:4f}s".format(str(ii+1)+"/"+str(N), time.time()-starttime))
        ii += 1
    print("Finished total computation in: {:4f}s".format(time.time()-totStartTime))
#    plt.xlabel("a1/a2")
#    plt.ylabel("c1-c2" )
#    plt.xlim([0,1])
#    plt.ylim([0,.1]) 
#    plt.show()


    # setup the information for saving histograms of good candidates
    fileroot = filename[8:-4]
    # create the folder for each cluster
    print("Creating folder for halo: {}".format(folder))
    if not os.path.exists("./halos/{}".format(folder)):
        print("Folder does not exist")
        os.makedirs("./halos/{}".format(folder))

    # make the fs array a 2d array
    fs = np.reshape(fs, (nPhi, nTheta))
    amplitudes = np.reshape(amplitudes, (nPhi, nTheta))
    widths = np.reshape(widths, (nPhi, nTheta))

    

    for theta in range(nTheta):
        for phi in range(nPhi):

            # create the figure
            histFig = plt.figure()
    
            f = fs[theta][phi]
            sigma = widths[theta][phi]
            A = amplitudes[theta][phi]     
            t = thetas[theta]
            p = phis[phi]
            # check if the sigma ratios are good
            if sigma > 0.3:
                # check if the amplitudes are comparable
                if A > 0.5:
                    # check if the f value is good
                    if f > 0.05:
                        #this should be a good one
                        print("Good candidate at theta: {} phi: {}".format(p*180/np.pi, t*180/np.pi))
                        
                        # create the figure
                        plt.figure(figsize=(10,5))
                        # recalculate the histograms from the good angles
                        vert = get_icovertices(0, [t], [p])[0]
                        ax = plt.subplot(1,2,1)
                        f, param1, param2 = compute_histogram(np.array(rotate_galaxies(relPositions, vert).T), 
                            np.array(rotate_galaxies(velocities, vert).T), masses, vert, randoms)

   


                        ax2 = plt.subplot(1,2,2)
                        # create a scatter plot of the cluster that has favorable qualities
                        scatter(np.array(rotate_galaxies(relPositions, vert).T), 
                                np.array(rotate_galaxies(velocities, vert).T), masses, vert, randoms)

                        plt.savefig("./halos/{}/{}_t{:.0f}p{:.0f}_s.png".format(folder, fileroot, p*180/np.pi, 
                                    t*180/np.pi))
                        plt.close()



            # close the plot so they don't get redrawn on top of eachother
            plt.close()

    plt.close()
    fig = plt.figure()
#    ax = fig.gca(projection='3d')
#    thetas, phis = np.meshgrid(thetas, phis)
#    fs = np.clip(fs, 0, 999)
#    surf = ax.plot_wireframe(thetas, phis, fs)
#    ax.set_zlim(0, 2)
#    ax.set_xlabel(r'$\theta$')
#    ax.set_ylabel(r'$\phi$')
    plt.pcolor(thetas*180/np.pi, phis*180/np.pi, fs, vmin=0)
    plt.xlabel(r'$\theta$', fontsize=16)
    plt.ylabel(r'$\phi$', fontsize=16)

    plt.colorbar().set_label("f", fontsize=16)



    fig1 = plt.figure()
    plt.pcolor(thetas*180/np.pi, phis*180/np.pi, widths, vmin=0)
    plt.xlabel(r'$\theta$', fontsize=16)
    plt.ylabel(r'\phi$', fontsize=16)

    plt.colorbar().set_label("Fractional $\sigma$", fontsize=16)
    
    fig2 = plt.figure()
    plt.pcolor(thetas*180/np.pi, phis*180/np.pi, amplitudes, vmin=0)
    plt.xlabel(r'$\theta$', fontsize=16)
    plt.ylabel(r'\phi$', fontsize=16)

    plt.colorbar().set_label("Fractionsl $A$", fontsize=16)

    # show the figures of fractional gaussian parameters and f 
#    plt.show()




    return relPositions, velocities, masses, randoms


# return the uniformly distributed vertices on an icosahedron
def get_icovertices(Nsubdiv, thetas, phis):
    # first create a regular icosahedron
    vertList = np.array([])

    g = (1+np.sqrt(5))/2.
    # a regular icosahedron is made from 3 orthogonal golden rectangles
    #  it will be added manually here
    vertList = np.append(vertList, [1, g, 0])
    vertList = np.append(vertList, [-1, g, 0])
    vertList = np.append(vertList, [1, -g, 0])
    vertList = np.append(vertList, [-1, -g, 0])

    vertList = np.append(vertList, [g, 0, 1])
    vertList = np.append(vertList, [-g, 0, 1])
    vertList = np.append(vertList, [g, 0, -1])
    vertList = np.append(vertList, [-g, 0, -1])

    vertList = np.append(vertList, [0, 1, g])
    vertList = np.append(vertList, [0, 1, -g])
    vertList = np.append(vertList, [0, -1, g])
    vertList = np.append(vertList, [0, -1, -g])

    vertList = np.reshape(vertList, (-1,3))

#    fig3d = plt.figure()
#    ax3d = fig3d.add_subplot(111, projection='3d')

    # plot it for testing purposes
#    for row in vertList:
#        ax3d.scatter(row[0], row[1], row[2])
#
#    plt.show()
    vertList = np.array([])
    for theta in thetas:
        for phi in phis:
            
            #vertList = np.append(vertList, [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)])
            vertList = np.append(vertList, [theta, phi])
            
    vertList = np.reshape(vertList, (-1, 2))
    return vertList


# this will rotate the positions based on the ico_vertices function then 
#  call the procedure to make this histogram
def rotate_galaxies(positions, vertex):
    # find the angles of rotation
    theta = vertex[0]
#    if np.isnan(vertex[0]/vertex[1]): 
#        theta = np.pi/2. 
#    else:
#        theta = atan(vertex[0]/vertex[1]) 
    Rz = np.matrix([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])

    phi = vertex[1]
#    vertex = Rz*np.matrix(vertex).T
#    if np.isnan(vertex[2]/vertex[0]):
#        phi = np.pi/2.
#    else:
#        phi = atan(vertex[2]/vertex[0])
    Ry = np.matrix([[cos(phi), 0, sin(phi)], [0, 1, 0], [-sin(phi), 0, cos(phi)]])

    # define the rotation matrices
    Rx = np.matrix([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])

    #multiply the position matrix by the rotation matrices
    positions = Rz*np.matrix(positions).T
    positions = Ry*np.matrix(positions)
     # plot it for testing purposes
#    fig3d = plt.figure()
#    ax3d = fig3d.add_subplot(111, projection='3d')
#    for row in positions.T:
#        print(row[0,0])
#        ax3d.scatter(row[0,0], row[0,1], row[0,2])

#    plt.show()
    return positions


# this will compute the 2d positional scatter plot of the data
def scatter(positions, velocities, masses, vert, randoms):
    theta = vert[0]*180/np.pi
    phi = vert[1]*180/np.pi




    # sort the lists based on mass
    zippedList = list(zip(masses, positions, velocities, randoms))
    zippedList = sorted(zippedList, key = lambda x: x[0])
    positions = np.array([])
    velocities = np.array([])
    masses = np.array([])
    randoms = np.array([])
    for entry in zippedList:
        positions = np.append(positions, entry[1])  
        masses = np.append(masses, entry[0])  
        velocities = np.append(velocities, entry[2])  
        randoms = np.append(randoms, entry[3])  
        # the lists are now sorted
        
    # reshape the arrays to be 3xwhatever
    positions = np.reshape(positions, (-1,3)) 
    velocities = np.reshape(positions, (-1,3))
    # define the number of laes and lbgs that will be present in the histogram
    nLAE = 100
    nLBG = 50
    # we will always choose the most massive halo to be an LBG
    nLBG -= 1


    #the cutoff minimum mass for LBGs
    LBGcutoff = 10**11.1
    # index of the cutoff
    LBGcutoffind = np.argmin(abs(masses-LBGcutoff))
    
    #occupation fraction of LAEs
    LAEoccup = .1
    
    # define the distance/redshift of the center of the cluster
    peakMidZ = 3.078
    peakMidMpc = cmToMpc*DC(peakMidZ)

    # pick out the arrays that will be used to find the LBGs and LAEs
    # LBG arrays
    LBGpos = np.array([])
    LBGvels = np.array([])
    LBGmass=np.array([])
    LBGrands = np.array([])
    # LAE arrays
    LAEpos = np.array([])
    LAEvels = np.array([])
    LAEmass=np.array([])
    LAErands = np.array([])

    # loop through the full lists and pick out the n with the lowest random numbers
    #  then assign the results 
    for ii in range(nLBG):
        # pick out the next LBG
        LBG = np.argmin(randoms[1:LBGcutoff])    
        LBGpos = np.append(LBGpos, positions[LBG])
        LBGvels = np.append(LBGvels, velocities[LBG])
        LBGmass=np.append(LBGmass, masses[LBG])
        LBGrands = np.append(LBGrands, randoms[LBG])
        # remove the taken LBG
        positions = np.delete(positions, LBG, 0)  
        masses = np.delete(masses, LBG, 0)  
        velocities = np.delete(velocities, LBG, 0)  
        randoms = np.delete(randoms, LBG, 0)  
       

        # there are a different number of elements above the cutoff now
        LBGcutoff -= 1
        

    # loop through the LAEs
    for jj in range(nLAE):
        LAE = np.argmin(randoms[1:]) 
        # pick out the next LAE
        LAEpos = np.append(LAEpos, positions[LAE])
        LAEvels = np.append(LAEvels, velocities[LAE])
        LAEmass=np.append(LAEmass, masses[LAE])
        LAErands = np.append(LAErands, randoms[LAE])

        # remove the taken LAE
        positions = np.delete(positions, LAE, 0)  
        masses = np.delete(masses, LAE, 0)  
        velocities = np.delete(velocities, LAE, 0)  
        randoms = np.delete(randoms, LAE, 0)  
 
    # reshape arrays
    LBGpos = np.reshape(LBGpos, (-1, 3))
    LBGvels = np.reshape(LBGvels, (-1,3))
    LAEpos = np.reshape(LAEpos, (-1,3))
    LAEvels = np.reshape(LAEvels, (-1, 3))

    # these are the final arrays with all of the LAE/LBG galaxies       
    galPos = np.reshape(np.append(LAEpos, LBGpos), (-1, 3))
    galVel = np.reshape(np.append(LAEvels, LBGvels), (-1, 3))




    # distance/redshift for the redshift center of the cluster, roughly
    peakMidZ = 3.078
    peakMidMpc = cmToMpc*DC(peakMidZ)

    #define colormap
    peakmin = 3.075
    peakmax = 3.082
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




    for row, vel in zip(galPos, galVel):
        x = row[0]
        # calculate the redshift here 
        z = find_redshift(peakMidMpc+x)
        # apply the galaxy velocity correction
        z += (vel[0]/3e5)*(1+z)

        size=50
        plt.scatter(row[1]/(7.78e-3)/60, row[2]/(7.78e-3)/60.,
        vmin = 3.07, vmax=3.083,
        c=z,
        cmap=colorMap, s=size, linewidth=.4)



    plt.xlabel("[arcmin]", fontsize=16)
    plt.ylabel("[arcmin]", fontsize=16)
    plt.colorbar().set_label("z", fontsize=16)





# given a set of data points compute the redshift histogram
def compute_histogram(positions, velocities, masses, vert, randoms):
    theta = vert[0]*180/np.pi
    phi = vert[1]*180/np.pi
    # we will only really care about the x-positions

    # sort the lists based on mass
    zippedList = list(zip(masses, positions, velocities, randoms))
    zippedList = sorted(zippedList, key = lambda x: x[0])
    positions = np.array([])
    velocities = np.array([])
    masses = np.array([])
    randoms = np.array([])
    for entry in zippedList:
        positions = np.append(positions, entry[1])  
        masses = np.append(masses, entry[0])  
        velocities = np.append(velocities, entry[2])  
        randoms = np.append(randoms, entry[3])  
        # the lists are now sorted
        
    # reshape the arrays to be 3xwhatever
    positions = np.reshape(positions, (-1,3)) 
    velocities = np.reshape(positions, (-1,3))
    # define the number of laes and lbgs that will be present in the histogram
    nLAE = 100
    nLBG = 50
    # we will always choose the most massive halo to be an LBG
    nLBG -= 1


    #the cutoff minimum mass for LBGs
    LBGcutoff = 10**11.1
    # index of the cutoff
    LBGcutoffind = np.argmin(abs(masses-LBGcutoff))
    
    #occupation fraction of LAEs
    LAEoccup = .1
    
    # define the distance/redshift of the center of the cluster
    peakMidZ = 3.078
    peakMidMpc = cmToMpc*DC(peakMidZ)

    # pick out the arrays that will be used to find the LBGs and LAEs
    # LBG arrays
    LBGpos = np.array([])
    LBGvels = np.array([])
    LBGmass=np.array([])
    LBGrands = np.array([])
    # LAE arrays
    LAEpos = np.array([])
    LAEvels = np.array([])
    LAEmass=np.array([])
    LAErands = np.array([])

    # loop through the full lists and pick out the n with the lowest random numbers
    #  then assign the results 
    for ii in range(nLBG):
        # pick out the next LBG
        LBG = np.argmin(randoms[1:LBGcutoff])    
        LBGpos = np.append(LBGpos, positions[LBG])
        LBGvels = np.append(LBGvels, velocities[LBG])
        LBGmass=np.append(LBGmass, masses[LBG])
        LBGrands = np.append(LBGrands, randoms[LBG])
        # remove the taken LBG
        positions = np.delete(positions, LBG, 0)  
        masses = np.delete(masses, LBG, 0)  
        velocities = np.delete(velocities, LBG, 0)  
        randoms = np.delete(randoms, LBG, 0)  
       

        # there are a different number of elements above the cutoff now
        LBGcutoff -= 1
        

    # loop through the LAEs
    for jj in range(nLAE):
        LAE = np.argmin(randoms[1:]) 
        # pick out the next LAE
        LAEpos = np.append(LAEpos, positions[LAE])
        LAEvels = np.append(LAEvels, velocities[LAE])
        LAEmass=np.append(LAEmass, masses[LAE])
        LAErands = np.append(LAErands, randoms[LAE])

        # remove the taken LAE
        positions = np.delete(positions, LAE, 0)  
        masses = np.delete(masses, LAE, 0)  
        velocities = np.delete(velocities, LAE, 0)  
        randoms = np.delete(randoms, LAE, 0)  
 
    # reshape arrays
    LBGpos = np.reshape(LBGpos, (-1, 3))
    LBGvels = np.reshape(LBGvels, (-1,3))
    LAEpos = np.reshape(LAEpos, (-1,3))
    LAEvels = np.reshape(LAEvels, (-1, 3))

    # these are the final arrays with all of the LAE/LBG galaxies       
    galPos = np.reshape(np.append(LAEpos, LBGpos), (-1, 3))
    galVel = np.reshape(np.append(LAEvels, LBGvels), (-1, 3))

    

    # array to hold the z data
    zs = np.array([])

    # loop through all rows in the position matrix and calculate the redshift
    #  of each object
    for row, vel in zip(galPos, galVel):

        x = row[0]
        # calculate the redshift here 
        z = find_redshift(peakMidMpc+x)
        # apply the galaxy velocity correction
        z += (vel[0]/3e5)*(1+z)

        zs = np.append(zs, z)
            
#    plt.savefig("./halo_images/"+str(time.time())+".png")
#    plt.close()
    fitParam,  chisqr = fit_histogram(zs, 2)
#    plt.scatter((min([a1,a2])/max([a1,a2])), abs(c1-c2), s=100/chisqr)
    if chisqr < 101:
        x = np.linspace(3.0,3.2, 1000)
        Gauss = twoGauss(x, *fitParam)
        #plt.plot(x, Gauss, linewidth=3)
#        plt.hist(zs, range=(3.07, 3.09), bins=20)
        plt.hist(zs, range=(3.05, 3.12), bins=30)
#        plt.ylim([0, 25])
#        plt.savefig("./halo_movie/"+str(theta)+"_"+str(phi)+".png")
#        plt.show()
#        plt.close()
    return fcalc(zs)






# this will calculate the redshift given a comoving distance
def find_redshift(d, OmegaM=0.308, OmegaL=0.692):
    starttime = time.time()
    OmegaK = 1-OmegaM-OmegaL
    dz = 0.001
    z = 3.05 
    dCalculated = 0
    while (dCalculated<d): 
        z += dz
        zarr = np.linspace(0,z, 200)
        
#        Ez = lambda z: 1./np.sqrt(OmegaM*(1+z)**3+OmegaK*(1+z)**2+OmegaL)
        Ez = 1./np.sqrt(OmegaM*(1+zarr)**3+OmegaK*(1+zarr)**2+OmegaL)

        
        dCalculated = cmToMpc*DH()*scipy.integrate.simps(Ez, zarr)

    return z



# this will fit the histogram to a double gaussian
def fit_histogram(zs, npeaks):
    
    #find histogram statistics/fits
    hist, bin_edges = np.histogram(zs, range=(3.05, 3.12), bins=30)
    #find the bin centers
    bin_centers = [0.5*(bin_edges[ii]+bin_edges[ii+1]) for ii in range(len(bin_edges)-1)]
    #find a fit to two gaussians
    try:
        if npeaks == 2:
            fitParam, fitCov = curve_fit(twoGauss, bin_centers, hist)
            a1 = fitParam[0]
            a2 = fitParam[3]
            c1 = fitParam[1]
            c2 = fitParam[4]
            chisqr = 0
        elif npeaks == 1:
            fitParam, fitCov = curve_fit(Gauss, bin_centers, hist)
            chisqr=0
            
        for ii in range(len(bin_centers)):
            if npeaks == 2:
                chisqr += (hist[ii]-twoGauss(bin_centers[ii], *fitParam))**2/twoGauss(bin_centers[ii], *fitParam)
            elif npeaks == 1:
                chisqr += (hist[ii]-Gauss(bin_centers[ii], *fitParam))**2/Gauss(bin_centers[ii], *fitParam)

    except RuntimeError:
        print("Error, no fit was found")
        fitParam = [-999, -999, -999, -999, -999, -999]
        chisqr = 100
    

    if chisqr>100:
        chisqr=100 
    return fitParam, chisqr




# function to calculate f
def fcalc(zs):
    

    # get the chisquared
    # chi1 corresponds to double peaked, chi2 corresponds to single peaked
    param2, chi2 = fit_histogram(zs, 1) 
    param1, chi1 = fit_histogram(zs, 2)

    # degrees of freedom in the fits
    dof2 = 3
    dof1 = 6
    f = (( chi1 - chi2 ) / (dof1 - dof2)) / (chi2 / dof2)

    if chi1 >= 100 or chi2 >= 100:
        f = 0

    zarr = np.linspace(3.06, 3.09, 100)
    plt.plot(zarr, twoGauss(zarr, *param1), 'k', linewidth=2)
    return f, param1, param2



# a function that will plot the halo mass function 
def halo_mass_function(masses, randoms):
    mp = 8.6e8
    #the cutoff minimum mass for LBGs
    LBGcutoff = 10**11.1
    
    #occupation fraction of LAEs
    LAEoccup = .1
#    plt.figure()
    # change the masses array to account for missing LAEs
    galMasses = np.array([])
    for mass, r in zip(masses, randoms):
        if mass < LBGcutoff:
            if r < LAEoccup:
                galMasses = np.append(galMasses, mass)
        else:
                galMasses = np.append(galMasses, mass)


    
    plt.hist(masses, bins=np.logspace(10, 13, 30), histtype="step", color='black', linestyle='dashed', label="Ntot={}".format(masses.size))
    plt.hist(galMasses, bins=np.logspace(10, 13, 30), histtype="step", color='black', label="Ngalaxies={}".format(galMasses.size))
    plt.ylim([1,1e3])
    plt.gca().set_xscale("log")
    plt.gca().set_yscale("log")
    plt.xlabel("Mass [M$_\odot$]", fontsize=16)
    plt.ylabel("N", fontsize=16)
    plt.legend(loc='upper right')







if __name__=="__main__":
    # this will split up each of the individual descendant halos, and willl result
    #  in the random numbers being re-generated
#    split_query("./halos/MillenniumSQL_Full.dat")
    
#    relPositions, velocities, masses, randoms = read_halos(0,"./halos/MillenniumSQL_Full_0.dat")

#    # loop through each of the massive halos
    for ii in range(18):
        relPositions, velocities, masses, randoms = read_halos(ii,"./halos/MillenniumSQL_Full_{}.dat".format(str(ii)))
#        halo_mass_function(masses, randoms)
#        # plot the LBG cutoff
#        plt.plot([10**11.1, 10**11.1], [0, 10**4], 'k--', linewidth=2)
#
#        plt.savefig("./HMFs/Halo_{}.png".format(ii), dpi=400)
#        plt.close()






    
#    for ii in range(18):
#        relPositions, velocities, masses, randoms = read_halos("./halos/MillenniumSQL_Full_{}.dat".format(str(ii)))
