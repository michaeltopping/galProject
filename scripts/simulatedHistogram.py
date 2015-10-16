import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from math import *
from cosmology import *
import scipy
import time
from mpl_toolkits.mplot3d import Axes3D



#Constants and conversions
cmToMpc = 3.241e-25
cmToGpc = 3.241e-28




#define two gaussians
def twoGauss(x, a1, c1, s1, a2, c2, s2):
    if s1 ==0 or s2 == 0:
        return 99999.9
    else:
        return a1*np.exp(-(x-c1)**2/(2*s1**2))+a2*np.exp(-(x-c2)**2/(2*s2**2))




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

        outFile.write("{}   {}    {}    {}    {}    {}    {}\n".format(row[0], row[18], row[19],
                                                    row[20], row[25], row[26], row[27]))
        
    


# read in the halo data from one file
def read_halos(filename):
    # first read in the data from one of the files
    data = np.genfromtxt(filename, dtype=str)


    positions = np.array([])
    velocities = np.array([])
    for row in data:
        # get the line in the positions array that will be added
        newline = [float(row[1]), float(row[2]), float(row[3])]
        newVelLine = [float(row[4]), float(row[5]), float(row[6])]
        positions = np.append(positions, newline)
        velocities = np.append(velocities, newVelLine)
        

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


    vertList = get_icovertices(0) 
    N = len(vertList)
    ii = 0
    totStartTime = time.time()
    for vert in vertList:
        starttime = time.time()
        compute_histogram(np.array(rotate_galaxies(relPositions, vert).T), 
                            np.array(rotate_galaxies(velocities, vert).T), vert)
        print("Finished raytrace: {} in {:4f}s".format(str(ii)+"/"+str(N), time.time()-starttime))
        ii += 1
    print("Finished total computation in: {:4f}s".format(time.time()-totStartTime))
    plt.xlabel("a1/a2")
    plt.ylabel("c1-c2" )
    plt.xlim([0,1])
    plt.ylim([0,.1]) 
    plt.show()


# return the uniformly distributed vertices on an icosahedron
def get_icovertices(Nsubdiv):
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
    theta = 0.
    phi = 0.
    nTheta = 1
    nPhi = 360
    for itheta in range(nTheta):
        for iphi in range(nPhi):
            
            #vertList = np.append(vertList, [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)])
            vertList = np.append(vertList, [theta, phi])
            phi += 2*np.pi/nPhi
            
        phi = 0
        theta += 2*np.pi/nTheta
    vertList = np.reshape(vertList, (-1, 2))
    print(vertList)
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


# given a set of data points compute the redshift histogram
def compute_histogram(positions, velocities, vert):
    theta = vert[0]*180/np.pi
    phi = vert[1]*180/np.pi
    # we will only really care about the x-positions
    
    # define the distance/redshift of the center of the cluster
    peakMidZ = 3.078
    peakMidMpc = cmToMpc*DC(peakMidZ)

    # array to hold the z data
    zs = np.array([])

    # loop through all rows in the position matrix and calculate the redshift
    #  of each object
    for row, vel in zip(positions, velocities):
        x = row[0]
        # calculate the redshift here 
        z = find_redshift(peakMidMpc+x)
        # apply the galaxy velocity correction
        z += (vel[0]/3e5)*(1+z)
        zs = np.append(zs, z)
#    plt.savefig("./halo_images/"+str(time.time())+".png")
#    plt.close()
    fitParam,  chisqr = fit_histogram(zs)
    print(chisqr)
#    plt.scatter((min([a1,a2])/max([a1,a2])), abs(c1-c2), s=100/chisqr)
    if chisqr < 100:
        x = np.linspace(3.0,3.2, 1000)
        Gauss = twoGauss(x, *fitParam)
        #plt.plot(x, Gauss, linewidth=3)
        plt.hist(zs, range=(3.07, 3.09), bins=20)
        plt.xlim([3.07, 3.09])
        plt.ylim([0,15])
        plt.savefig("./halo_movie/"+str(theta)+"_"+str(phi)+".png")
        plt.close()


# this will calculate the redshift given a comoving distance
def find_redshift(d, OmegaM=0.308, OmegaL=0.692):
    starttime = time.time()
    OmegaK = 1-OmegaM-OmegaL
    dz = 0.001
    z = 2.9 
    dCalculated = 0
    while (dCalculated<d): 
        z += dz
        zarr = np.linspace(0,z, 200)
        
#        Ez = lambda z: 1./np.sqrt(OmegaM*(1+z)**3+OmegaK*(1+z)**2+OmegaL)
        Ez = 1./np.sqrt(OmegaM*(1+zarr)**3+OmegaK*(1+zarr)**2+OmegaL)

        
        dCalculated = cmToMpc*DH()*scipy.integrate.simps(Ez, zarr)

    return z



# this will fit the histogram to a double gaussian
def fit_histogram(zs):
    
    #find histogram statistics/fits
    hist, bin_edges = np.histogram(zs, range=(3.07, 3.09), bins=20)
    #find the bin centers
    bin_centers = [0.5*(bin_edges[ii]+bin_edges[ii+1]) for ii in range(len(bin_edges)-1)]
    #find a fit to two gaussians
    try:
        fitParam, fitCov = curve_fit(twoGauss, bin_centers, hist)
        a1 = fitParam[0]
        a2 = fitParam[3]
        c1 = fitParam[1]
        c2 = fitParam[4]
        chisqr = 0
        for ii in range(len(bin_centers)):
            chisqr += (hist[ii]-twoGauss(bin_centers[ii], *fitParam))**2/twoGauss(bin_centers[ii], *fitParam)

    except RuntimeError:
        print("Error, no fit was found")
        fitParam = [-999, -999, -999, -999, -999, -999]
        chisqr = 100
    

    if chisqr>100:
        chisqr=100 
    return fitParam, chisqr






if __name__=="__main__":
    split_query("./halos/testSQL.dat")
    read_halos("./halos/testSQL_0.dat")
