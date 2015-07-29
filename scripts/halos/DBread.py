import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#create a class for a cluster, that will hold many halos
class cluster():
    def __init__(self):
        self.halos = []



#create a procedure that returns the average of a field in an array of dicts
def dictAvg(arr, field):
    n = 0
    tot = 0
    for d in arr:
        tot += d[field]  
        n += 1  
    return tot/n 

#create a procedure that returns the standard deviation of a field in an 
# array of dicts
def dictStd(arr, field):
    n = 0
    tot = 0
    avg = dictAvg(arr, field)
    for d in arr:
        tot += (d[field]-avg)**2

    return np.sqrt(tot/(n-1))




#Read in the result from a millenium SQL query and sort through the data
def DBread(filename):
    #Read in the data table 
    data = np.genfromtxt(filename, dtype=str, skip_header=55, delimiter=',')

    #create the array that will hold split halo information
    clusters = {}
    #define the first descendent halo
    currentDesHalo = data[0][0]
    clusters[currentDesHalo] = cluster()
    #loop through each row
    for row in data:
        #check if this cluster is in the dictionary
        if not row[0] in clusters:
            #if the cluster does not exist, create it
            clusters[row[0]] = cluster()
        #add the halo to the cluster
        clusters[row[0]].halos.append({'x':float(row[18]), 
          'y':float(row[19]), 'z':float(row[20])})
    return clusters

#plot one of the clusters
def plotCluster(cluster, marker='o'):
    fig3d = plt.figure()
    ax3d = fig3d.add_subplot(111, projection='3d')
    avgx = dictAvg(cluster.halos, 'x')
    avgy = dictAvg(cluster.halos, 'y')
    avgz = dictAvg(cluster.halos, 'z')
    for halo in cluster.halos:
        ax3d.scatter(halo['x']-avgx, halo['y']-avgy, halo['z']-avgz, marker=marker)


#main
if __name__ == "__main__":
    starttime = time.time()
    print("Running as MAIN")
    clusters = DBread("testSQL.dat")
    clusterNames = ['1000019001111', '2000020000000', '3000042000000']
    plotCluster(clusters[clusterNames[0]], marker='v')
    plotCluster(clusters[clusterNames[1]], marker='o')
    
    print("Finished in: {:.3f} seconds".format(time.time()-starttime))
    plt.show()    
