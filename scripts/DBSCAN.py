#implementation of the DBSCAN Algorithm

import numpy as np
import matplotlib.pyplot as plt
import random
from math import *


#create some classes to be used to hold data
class Cluster:
    
    def __init__(self):
        self.points = []


class point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.visited = 0
        self.noise = 0
        self.inCluster = 0

    def regionQuery(self, D, eps):
        #return all points within eps
        inRegion = []
        for P in D:
            distance = sqrt((P.x-self.x)**2+(P.y-self.y)**2 +(P.z-self.z)**2)  
            if distance <= eps:
                inRegion.append(P)
        return inRegion


#define the main function DBSCAN
#takes input: D - list of points
#             eps - distance between connected points
#             minPts - minimum number of points to be considered a core point
def DBSCAN(D, eps, minPts):
    #create a set of points for classes
    pts = []
    noisePts = []
    for row in D:
        pts.append(point(row[0], row[1], row[2]))

    clusters = []
    #loop through each point in the list
    for P in pts:
        #check if this point has been visited
        if P.visited:
            continue
        #this point has now been visited
        P.visited=1
        #get a list of all points within eps of the point
        neighborPts = P.regionQuery(pts, eps)
        #check if this is noise
        if len(neighborPts) < minPts:
            P.noise = 1
            noisePts.append(P)
        else:
            clusters.append(Cluster())
            cluster = clusters[-1]
            expandCluster(P, neighborPts, cluster, eps, minPts, pts)
    
    return clusters, noisePts


#expand the cluster
def expandCluster(P, neighborPts, cluster, eps, minPts, pts):
    P.inCluster = 1
    cluster.points.append(P)
    for Q in neighborPts:
        if not Q.visited:
            Q.visited = 1
            neighborQts = Q.regionQuery(pts, eps)
            if len(neighborQts) >= minPts:
                neighborPts +=neighborQts
        if not Q.inCluster:
            Q.inCluster = 1
            cluster.points.append(Q)




#expandCluster will



if __name__ == "__main__":
    random.seed(6563)
    #create a list of random gaussian points
    group1 = np.array([])
    for ii in range(1000):
        #create random x and y coords with mean 0 and sigma 1
        x = random.gauss(0, 1)
        y = random.gauss(0, 1)
        group1 = np.append(group1, [x, y])
        x = random.gauss(8, 1)
        y = random.gauss(7, 1)
        group1 = np.append(group1, [x, y])
    #add in random points
    for jj in range(100):
        group1 = np.append(group1, [10*random.random(), 10*random.random()])
    group1 = group1.reshape((-1, 2))


    clusters, noise = DBSCAN(group1, 0.2, 3)

    colors = ['k', 'b', 'r', 'y', 'c']
    #loop through and plot the points
    for ii in range(len(clusters)):
        for P in clusters[ii].points:
            plt.plot(P.x, P.y, '{}^'.format(colors[ii%4]))

    #plot the noise points
    for P in noise:
        plt.plot(P.x, P.y, 'k.')

    plt.show()

