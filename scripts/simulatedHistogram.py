import numpy as np
import matplotlib.pyplot as plt
import os
from math import *
from cosmology import *
import scipy
import time
from mpl_toolkits.mplot3d import Axes3D



#Constants and conversions
cmToMpc = 3.241e-25
cmToGpc = 3.241e-28



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

        outFile.write("{}   {}    {}    {}\n".format(row[0], row[18], row[19],
                                                    row[20]))
        
    


# read in the halo data from one file
def read_halos(filename):
    # first read in the data from one of the files
    data = np.genfromtxt(filename, dtype=str)


    positions = np.array([])
    for row in data:
        # get the line in the positions array that will be added
        newline = [float(row[1]), float(row[2]), float(row[3])]
        positions = np.append(positions, newline)

    # reshape the array to be a 3xN array
    positions = np.reshape(positions, (-1,3))
    
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
    for vert in vertList:
        compute_histogram(np.array(rotate_galaxies(relPositions, vert).T))

    



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
    for itheta in range(10):
        for iphi in range(10):
            
            #vertList = np.append(vertList, [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)])
            vertList = np.append(vertList, [theta, phi])
            phi += 2*np.pi/10.
            
        theta += 2*np.pi/10.
    vertList = np.reshape(vertList, (-1, 2))
    return vertList


# this will rotate the positions based on the ico_vertices function then 
#  call the procedure to make this histogram
def rotate_galaxies(positions, vertex):
    print("Vertex: {}".format(vertex))
    # find the angles of rotation
    theta = vertex[0]
#    if np.isnan(vertex[0]/vertex[1]): 
#        theta = np.pi/2. 
#    else:
#        theta = atan(vertex[0]/vertex[1]) 
    print("Theta: {}".format(theta))

    Rz = np.matrix([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])

    phi = vertex[1]
#    vertex = Rz*np.matrix(vertex).T
#    if np.isnan(vertex[2]/vertex[0]):
#        phi = np.pi/2.
#    else:
#        phi = atan(vertex[2]/vertex[0])
    print("phi: {}".format(phi))
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
    print("--------------")
    print(positions.T)
    return positions


# given a set of data points compute the redshift histogram
def compute_histogram(positions):
    # we will only really care about the x-positions
    
    # define the distance/redshift of the center of the cluster
    peakMidZ = 3.078
    peakMidMpc = cmToMpc*DC(peakMidZ)

    # array to hold the z data
    zs = np.array([])

    # loop through all rows in the position matrix and calculate the redshift
    #  of each object
    for row in positions:
        x = row[0]
        # calculate the redshift here 
        z = find_redshift(peakMidMpc+x)
        zs = np.append(zs, z)
    plt.hist(zs)
    plt.savefig("./halo_images/"+str(time.time())+".png")
    plt.close()


# this will calculate the redshift given a comoving distance
def find_redshift(d, OmegaM=0.308, OmegaL=0.692):
    starttime = time.time()
    OmegaK = 1-OmegaM-OmegaL
    dz = 0.0005
    z = 2.5 
    dCalculated = 0
    while (dCalculated<d): 
        z += dz
        zarr = np.linspace(0,z, 200)
        
#        Ez = lambda z: 1./np.sqrt(OmegaM*(1+z)**3+OmegaK*(1+z)**2+OmegaL)
        Ez = 1./np.sqrt(OmegaM*(1+zarr)**3+OmegaK*(1+zarr)**2+OmegaL)

        
        dCalculated = cmToMpc*DH()*scipy.integrate.simps(Ez, zarr)

    return z


if __name__=="__main__":
    read_halos("./halos/testSQL_0.dat")
