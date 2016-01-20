import numpy as np
import os
from math import *
import matplotlib.pyplot as plt
from matplotlib import rc
from treeClasses import *
from cosmology import *
import random





# change the fonts on the plots

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# physical constants
c = 3e10
sToYr = (1./3.15e7)



# define a function to break up the query into the different cluster trees
def tree_break(filename):
    data = np.genfromtxt(filename, dtype=str, skip_header=57, delimiter=',')
    

    # create a dictionary to map the descendent halo to a cluster/tree
    desTocluster = {}
    nCluster = 0

    # create the list of files
    files = []

    # loop through each row of the data
    for row in data:
    
        # check if the current cluster is in the dictionary
        if not row[0] in desTocluster:
            # add the cluster to the dictionary
            desTocluster[row[0]] = nCluster

            # open the file for appending
            files.append(open('./halos/tree_{}.dat'.format(nCluster), 'a'))    
            print("Created file for tree/cluster: {}".format(nCluster))

            # move on to the next cluster for the future
            nCluster += 1

        fileNum = row[0]

        # format the string for writing to file
        charReplace = "'[]\n"
        row = np.array_str(row)
        for c in charReplace:
            row = row.replace(c, "")
           
        # add the row to the appropriate file
        files[desTocluster[fileNum]].write(row+"\n")




# create a function to break up a tree into snapshots
def snapnum_break(filename):

    data = np.genfromtxt(filename, dtype=str)


    # break up the files into redshifts
    # create a folder for the tree
    if not os.path.exists("./halos/treeZs_0"):
        os.makedirs("./halos/treeZs_0")

    files = []
    # create a file for each snapnum
    for ii in range(64): 
       files.append(open("./halos/treeZs_0/{}.dat".format(ii), 'a')) 


    for row in data:
        snapnum = int(row[13])
        # format the string for writing to file
        charReplace = "'[]\n"
        row = np.array_str(row)
        for c in charReplace:
            row = row.replace(c, "")


        files[snapnum].write(row+"\n")
 



# given a filename, make a halo tree
def make_tree(filename):

    data = np.genfromtxt(filename, dtype=str)



    # create the dictionary of halos
    halos = {}

    for row in data:

        # check if the halo is in the dictionary
        if not row[2] in halos:
            halos[row[2]] = halo(row) 
        # if the halo is in the dictionary
        else:
            print("Error: halo {} already exists in the dicionary".format(row[2]))


        # organize the halos so that each halo points to a descendent
#        for h in halos:
#            hObj = halos[h]
#           
#            if int(hObj.desid) > 0 and (hObj.desid in halos):
#                hObj.set_des(halos[hObj.desid])
#            else:
#                hObj.set_des(0)


    snapnums = []
       
    # loop through each of the snapshots
    for ii in range(64):
        firstIds = []
        # read in the file for the particular snapshop
        data = np.genfromtxt("./halos/treeZs_0/{}.dat".format(ii))

        # create the dictionary of halos
        snapHalos = {}

        for row in data:

            # check if the halo is in the dictionary
            if not row[2] in snapHalos:
                snapHalos[row[2]] = halo(row) 


        snapnums.append(snapnum(ii)) 
        # looo pthrough each of the halos in this snapshot
        for row in data:

            firstFofId = int(row[10])
            if not firstFofId in firstIds:
                firstIds.append(firstFofId)
                snapnums[ii].create_fof(snapHalos, firstFofId)

            
        
        print("Finished shapshot {} with {} fof_groups, containing {} halos".format(ii, len(snapnums[ii].fofGroups), np.shape(data)[0]))


        
    # find the maximum number of halos:
    maxhalos = 0
    maxgroups = 0
    for snap in snapnums:
        snap.remove_small_groups(3*10**12)
        snap.get_nHalos()
        if snap.nHalos > maxhalos:
            maxhalos = snap.nHalos            
        if snap.nGroups > maxgroups:
            maxgroups = snap.nGroups

       

    # arrange the xs for each snapnum
    for snap in snapnums:
        # find the number of halos in the snapshot
        snap.arrange_xs(maxgroups)


    for ii in range(len(snapnums)):
        snapnums[ii].plot_snapshot()
        if ii < 63: 
            snapnums[ii].plot_prog(snapnums[ii+1])
 
    plt.ylim([0, 3e10])
    # set the axes labels
    ax1 = plt.gca()
    ax1.set_ylabel('Lookback Time')
    ax2 = ax1.twinx()
     
    ax2Zs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    ax2Ts = [DC(z)*sToYr/c for z in ax2Zs]
    ax2.set_yticks(ax2Ts)
    ax2.set_yticklabels(ax2Zs)
    ax2.set_ylabel('z')
    plt.show()


#    # loop through the file until all of the halos are empty
#    halokeys = list(halos.keys())
#    while halos:
#        for key in halokeys:
#            if key in halos:
#                print("Create fof group number {}".format(ifof))
#                fof = fof_group(ifof)
#                ifof += 1
#                halos = fof.fill_halos(halos, halos[key].firstProgId)
#                if key in halos:
#                    del halos[key]
#                print(len(halos))
#                print(key)
#                print(len(halokeys))
#                if len(fof.halos) > 1:
#                    fofs.append(fof)
#        halokeys = list(halos.keys())
#        
#
#    print("Plotting Tree")
#    for f in fofs:
#        mass = 0
#        for h in f.halos:
#            mass += h.m
#            print(h.m) 
#        z = f.halos[0].z
#        plt.plot(.5, DC(z)*sToYr/c, 'o', markersize=(max(sqrt(mass), 5e6))*5e-7)
#    plt.show()


    # go through and plot all of the halos
#    for h in halos:
#        hObj = halos[h]
#        plt.plot(random.random(), DC(hObj.z)*sToYr/c, 'o')
#
#    plt.show()

        



    







# main
if __name__ == "__main__":
#    tree_break('./halos/MR_Halos.dat')
#    snapnum_break('./halos/tree_0.dat')
    make_tree("./halos/tree_0.dat")
