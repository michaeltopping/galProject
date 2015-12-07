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

from matplotlib import rc


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


#Constants and conversions
cmToMpc = 3.241e-25
cmToGpc = 3.241e-28







# read in the halo data from one file
def read_halos(filename):
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
    
    return relPositions






# plot the visualization in 3d
def plot3D(positions, angle):
    # create the figure
    fig3d = plt.figure()
    ax3d = fig3d.add_subplot(111, projection='3d')


    for pos in positions:
        ax3d.scatter(pos[0], pos[1], pos[2])

    plt.show()

# this program will create a 3d visualization of the halos from the Millennium simulation
if __name__=="__main__":
    
    plot3D(read_halos("./halos/MillenniumSQL_Full_0.dat"), 0)



