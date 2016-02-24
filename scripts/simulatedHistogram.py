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
import pickle


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




# find all of the fof_groups in the snapshot
def find_fof(filename):

    # read in the result of the sql query
    data = np.genfromtxt(filename, dtype=str, skip_header=55, delimiter=',')

    
    # the list will be a dictionary with fofId:[haloids]
    groups = {}
        

    fofIDs = {}
    nfof = 0
    currentDesId = data[0][0]
    nCluster = 0
    # loop through the rows to create the list of first halos in fof groups
    for row in data:
        print(nCluster, len(groups)) 
        if currentDesId == row[0]:
            if not row[10] in fofIDs:
                fofIDs[row[10]] = nfof
                nfof += 1
        else:
            # loop through all of the rows
            for entry in data:
                if entry[0] == currentDesId:
                    groups[entry[2]] = fofIDs[entry[10]]


            pickle.dump((fofIDs, groups), open("./halos/fof_groupIDs_{}.p".format(nCluster), "wb"))
            groups = {}
            nfof = 0
            fofIDs = {}
            
            currentDesId = row[0]
            fofIDs[row[10]] = nfof
            nfof += 1
            nCluster += 1



    for entry in data:
        if entry[0] == currentDesId:
            groups[entry[2]] = fofIDs[entry[10]]

      
        
    pickle.dump((fofIDs,groups), open("./halos/fof_groupIDs_{}.p".format(nCluster), "wb"))
        
        









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

        outFile.write("{}   {}    {}    {}    {}    {}    {}    {}    {}    {}\n".format(row[0], row[19], row[20],row[21], row[23], row[24], row[25], row[15], random.random(), row[2]))
        
    


# read in the halo data from one file
def read_halos(folder, filename):
    # first read in the data from one of the files
    data = np.genfromtxt(filename, dtype=str)
    mp = 8.6e8


    positions = np.array([])
    velocities = np.array([])
    masses = np.array([])
    randoms = np.array([])
    haloIDs = np.array([])
    for row in data:
        # get the line in the positions array that will be added
        newline = [float(row[1]), float(row[2]), float(row[3])]
        newVelLine = [float(row[4]), float(row[5]), float(row[6])]
        positions = np.append(positions, newline)
        velocities = np.append(velocities, newVelLine)
        masses = np.append(masses, float(row[7])*mp)
        randoms = np.append(randoms, float(row[8]))
        haloIDs = np.append(haloIDs, float(row[9]))
        

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


    nTheta = 10
    nPhi = 10
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

        # here we will find the galaxies that are within the most optimal field of veiw in the 
        #  millennium simulation
#        fovIndices = get_fov(np.array(rotate_galaxies(relPositions, vert).T))
        

        # note that the parameters for single and double peaked fits
        #  are switched here because I was to lazy to switch all following references
        #  to them. :(
        f, param2, param1 = compute_histogram(np.array(rotate_galaxies(relPositions, vert).T), 
                            np.array(rotate_galaxies(velocities, vert).T), masses, vert, randoms, haloIDs, folder)

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
#            if sigma > 0.3:
#                # check if the amplitudes are comparable
#                if A > 0.5:
#                    # check if the f value is good
#                    if f > 0.05:
                        #this should be a good one
            print("Good candidate at theta: {} phi: {}".format(p*180/np.pi, t*180/np.pi))
            
            # create the figure
            plt.figure(figsize=(10,10))
            # recalculate the histograms from the good angles
            vert = get_icovertices(0, [t], [p])[0]

            # here we will find the galaxies that are within the most optimal field of veiw in the 
            #  millennium simulation
            fovIndices = get_fov(np.array(rotate_galaxies(relPositions, vert).T))
            print("fovIndices: {}".format(fovindices))
 
            relPositions = relPositions[fovIndices]
            velocities = velocities[fovIndices]
            masses = masses[fovIndices]
            randoms = randoms[fovIndices]
            haloIDs = haloIDs[fovIndices]


            ax = plt.subplot(2,2,1)
            f, param1, param2 = compute_histogram(np.array(rotate_galaxies(relPositions, vert).T), 
                np.array(rotate_galaxies(velocities, vert).T), masses, vert, randoms, haloIDs, folder)




            ax2 = plt.subplot(2,2,2)
            # create a scatter plot of the cluster that has favorable qualities
            scatter(np.array(rotate_galaxies(relPositions, vert).T), 
                    np.array(rotate_galaxies(velocities, vert).T), masses, vert, randoms, haloIDs, folder)

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
def scatter(positions, velocities, masses, vert, randoms, haloIDs, ii):


    theta = vert[0]*180/np.pi
    phi = vert[1]*180/np.pi


    # get the list of fofIds for each halo
    fofIds, groups = pickle.load(open("./halos/fof_groupIDs_{}.p".format(ii), "rb"))
    print("Finished unpickling group data")


    # find the most massive fof groups
    fofMasses = {}
    # loop through all halos
    for groupId in groups:
        # find the mass of that halo
        
        mass = masses[np.where(haloIDs == int(groupId))]
        # add that mass onto the mass of the fof group
        if not groups[groupId] in fofMasses:
            fofMasses[groups[groupId]] = mass
        else:
            fofMasses[groups[groupId]] += mass

    

    nmaxs = 2
    maxs = [x[0] for x in sorted(fofMasses.items(), key = operator.itemgetter(1), reverse=True)[0:nmaxs]]
    print(maxs)
    print("Found maxs with mass:{},{}; {},{}".format(maxs[0], fofMasses[maxs[0]],maxs[1], fofMasses[maxs[1]]))


    
    


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
    nLAE = 116
    nLBG = 86


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
    LBGIDs = np.array([])
    # LAE arrays
    LAEpos = np.array([])
    LAEvels = np.array([])
    LAEmass=np.array([])
    LAErands = np.array([])
    LAEIDs = np.array([])


    # loop through all of the halos and add all of the ones that are in the 
    #  first or second most massive halos
#    deletes = np.array([])
#    for ii in range(randoms.size):
#        if ( groups[str(int(haloIDs[ii]))] == maxs[0] ):
#            LBG = ii
#            LBGpos = np.append(LBGpos, positions[LBG])
#            LBGvels = np.append(LBGvels, velocities[LBG])
#            LBGmass=np.append(LBGmass, masses[LBG])
#            LBGrands = np.append(LBGrands, randoms[LBG])
#            LBGIDs = np.append(LBGIDs, haloIDs[LBG])
#            deletes = np.append(deletes, LBG)
#            nLBG -= 1
#        
#        if ( groups[str(int(haloIDs[ii]))] == maxs[1] ):
#            LBG = ii
#            LBGpos = np.append(LBGpos, positions[LBG])
#            LBGvels = np.append(LBGvels, velocities[LBG])
#            LBGmass=np.append(LBGmass, masses[LBG])
#            LBGrands = np.append(LBGrands, randoms[LBG])
#            LBGIDs = np.append(LBGIDs, haloIDs[LBG])
#            deletes = np.append(deletes, LBG)
#            nLBG -= 1
#
#
#
#    # delete the LBGs that were used in the massive halos
#    for d in np.sort(deletes)[::-1]:
#        positions = np.delete(positions, d, 0)  
#        masses = np.delete(masses, d, 0)  
#        velocities = np.delete(velocities, d, 0)  
#        randoms = np.delete(randoms, d, 0)  
#        haloIDs = np.delete(haloIDs, d, 0)
 



    # loop through the full lists and pick out the n with the lowest random numbers
    #  then assign the results 
    for ii in range(nLBG):
        # check if there are enough left
        if randoms.shape[0] > 0:
            # pick out the next LBG
            LBG = np.argmin(randoms[1:LBGcutoff])    
            LBGpos = np.append(LBGpos, positions[LBG])
            LBGvels = np.append(LBGvels, velocities[LBG])
            LBGmass=np.append(LBGmass, masses[LBG])
            LBGrands = np.append(LBGrands, randoms[LBG])
            LBGIDs = np.append(LBGIDs, haloIDs[LBG])
            # remove the taken LBG
            positions = np.delete(positions, LBG, 0)  
            masses = np.delete(masses, LBG, 0)  
            velocities = np.delete(velocities, LBG, 0)  
            randoms = np.delete(randoms, LBG, 0)  
            haloIDs = np.delete(haloIDs, LBG, 0)
           

            # there are a different number of elements above the cutoff now
            LBGcutoff -= 1
        

    # loop through the LAEs
    for jj in range(nLAE):

        # check if there are enough galaxies left
        if randoms.shape[0] > 0:

            LAE = np.argmin(randoms[1:]) 
            # pick out the next LAE
            LAEpos = np.append(LAEpos, positions[LAE])
            LAEvels = np.append(LAEvels, velocities[LAE])
            LAEmass=np.append(LAEmass, masses[LAE])
            LAErands = np.append(LAErands, randoms[LAE])
            LAEIDs = np.append(LAEIDs, haloIDs[LAE])

            # remove the taken LAE
            positions = np.delete(positions, LAE, 0)  
            masses = np.delete(masses, LAE, 0)  
            velocities = np.delete(velocities, LAE, 0)  
            randoms = np.delete(randoms, LAE, 0)  
            haloIDs = np.delete(haloIDs, LAE, 0)
 
    # reshape arrays
    LBGpos = np.reshape(LBGpos, (-1, 3))
    LBGvels = np.reshape(LBGvels, (-1,3))
    LAEpos = np.reshape(LAEpos, (-1,3))
    LAEvels = np.reshape(LAEvels, (-1, 3))

    # these are the final arrays with all of the LAE/LBG galaxies       
    galPos = np.reshape(np.append(LAEpos, LBGpos), (-1, 3))
    galVel = np.reshape(np.append(LAEvels, LBGvels), (-1, 3))
    galIDs = np.append(LAEIDs, LBGIDs)




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




    for row, vel, haloId in zip(galPos, galVel, galIDs):
#        haloId = str(int(haloId)) 
#        # find the plotting parameters based on the fof group   
#        if groups[haloId] == maxs[0]:
#            sym = 's' 
#            size=70
#        elif groups[haloId] == maxs[1]:
#            sym = '^' 
#            size=70
#        else:
#            sym = 'o'
#            size=25
            

        x = row[0]
        # calculate the redshift here 
        z = find_redshift(peakMidMpc+x)
        # apply the galaxy velocity correction
        z += (vel[0]/3e5)*(1+z)

        plt.scatter(row[1]/(7.78e-3)/60, row[2]/(7.78e-3)/60., 
        vmin = 3.07, vmax=3.083,
        c=z, marker='o',
        cmap=colorMap, s=40, linewidth=.4)



    plt.xlabel("[arcmin]", fontsize=16)
    plt.ylabel("[arcmin]", fontsize=16)
    plt.colorbar().set_label("z", fontsize=16)





# given a set of data points compute the redshift histogram
def compute_histogram(positions, velocities, masses, vert, randoms, haloIDs, ii):



    theta = vert[0]*180/np.pi
    phi = vert[1]*180/np.pi


    # get the list of fofIds for each halo
    fofIds, groups = pickle.load(open("./halos/fof_groupIDs_{}.p".format(ii), "rb"))
    print("Finished unpickling group data")


    # find the most massive fof groups
    fofMasses = {}
    # loop through all halos
    for groupId in groups:
        # find the mass of that halo
        
        mass = masses[np.where(haloIDs == int(groupId))]
        # add that mass onto the mass of the fof group
        if not groups[groupId] in fofMasses:
            fofMasses[groups[groupId]] = mass
        else:
            fofMasses[groups[groupId]] += mass

     

    nmaxs = 2
    maxs = [x[0] for x in sorted(fofMasses.items(), key = operator.itemgetter(1), reverse=True)[0:nmaxs]]
    print("Maxs:{}".format(maxs))
    print("Found maxs with mass:{},{}; {},{}".format(maxs[0], fofMasses[maxs[0]],maxs[1], fofMasses[maxs[1]]))


    


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
    nLAE = 116
    nLBG = 86


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
    LBGIDs = np.array([])
    # LAE arrays
    LAEpos = np.array([])
    LAEvels = np.array([])
    LAEmass=np.array([])
    LAErands = np.array([])
    LAEIDs = np.array([])



    # loop through all of the halos and add all of the ones that are in the 
    #  first or second most massive halos
#    deletes = np.array([])
#    for ii in range(randoms.size):
#        if ( groups[str(int(haloIDs[ii]))] == maxs[0] ):
#            LBG = ii
#            LBGpos = np.append(LBGpos, positions[LBG])
#            LBGvels = np.append(LBGvels, velocities[LBG])
#            LBGmass=np.append(LBGmass, masses[LBG])
#            LBGrands = np.append(LBGrands, randoms[LBG])
#            LBGIDs = np.append(LBGIDs, haloIDs[LBG])
#            deletes = np.append(deletes, LBG)
#            nLBG -= 1
#        
#        if ( groups[str(int(haloIDs[ii]))] == maxs[1] ):
#            LBG = ii
#            LBGpos = np.append(LBGpos, positions[LBG])
#            LBGvels = np.append(LBGvels, velocities[LBG])
#            LBGmass=np.append(LBGmass, masses[LBG])
#            LBGrands = np.append(LBGrands, randoms[LBG])
#            LBGIDs = np.append(LBGIDs, haloIDs[LBG])
#            deletes = np.append(deletes, LBG)
#            nLBG -= 1
#
#
#
#    # delete the LBGs that were used in the massive halos
#    for d in np.sort(deletes)[::-1]:
#        positions = np.delete(positions, d, 0)  
#        masses = np.delete(masses, d, 0)  
#        velocities = np.delete(velocities, d, 0)  
#        randoms = np.delete(randoms, d, 0)  
#        haloIDs = np.delete(haloIDs, d, 0)
        
        

            
        




    # loop through the full lists and pick out the n with the lowest random numbers
    #  then assign the results 
    for ii in range(nLBG):
        # check if there are enough galaxies left
        if randoms.shape[0] > 0:
            # pick out the next LBG
            LBG = np.argmin(randoms[1:LBGcutoff])    
            LBGpos = np.append(LBGpos, positions[LBG])
            LBGvels = np.append(LBGvels, velocities[LBG])
            LBGmass=np.append(LBGmass, masses[LBG])
            LBGrands = np.append(LBGrands, randoms[LBG])
            LBGIDs = np.append(LBGIDs, haloIDs[LBG])
            # remove the taken LBG
            positions = np.delete(positions, LBG, 0)  
            masses = np.delete(masses, LBG, 0)  
            velocities = np.delete(velocities, LBG, 0)  
            randoms = np.delete(randoms, LBG, 0)  
            haloIDs = np.delete(haloIDs, LBG, 0)
           

            # there are a different number of elements above the cutoff now
            LBGcutoff -= 1
        

    # loop through the LAEs
    for jj in range(nLAE):
        # check if there are enough galaxies left
        if randoms.shape[0] > 0:
            LAE = np.argmin(randoms[1:]) 
            # pick out the next LAE
            LAEpos = np.append(LAEpos, positions[LAE])
            LAEvels = np.append(LAEvels, velocities[LAE])
            LAEmass=np.append(LAEmass, masses[LAE])
            LAErands = np.append(LAErands, randoms[LAE])
            LAEIDs = np.append(LAEIDs, haloIDs[LAE])

            # remove the taken LAE
            positions = np.delete(positions, LAE, 0)  
            masses = np.delete(masses, LAE, 0)  
            velocities = np.delete(velocities, LAE, 0)  
            randoms = np.delete(randoms, LAE, 0)  
            haloIDs = np.delete(haloIDs, LAE, 0)
     
    # reshape arrays
    LBGpos = np.reshape(LBGpos, (-1, 3))
    LBGvels = np.reshape(LBGvels, (-1,3))
    LAEpos = np.reshape(LAEpos, (-1,3))
    LAEvels = np.reshape(LAEvels, (-1, 3))




    # these are the final arrays with all of the LAE/LBG galaxies       
    galPos = np.reshape(np.append(LAEpos, LBGpos), (-1, 3))
    galVel = np.reshape(np.append(LAEvels, LBGvels), (-1, 3))
    galIDs = np.append(LAEIDs, LBGIDs)




    # array to hold the z data
    # also hold the arrays for the first and second most massive cluster
    firstMassive = np.array([])
    secondMassive = np.array([])
    lessMassive = np.array([])
    zs = np.array([])

    # loop through all rows in the position matrix and calculate the redshift
    #  of each object
    for row, vel, haloId in zip(galPos, galVel, galIDs):

        x = row[0]
        # calculate the redshift here 
        z = find_redshift(peakMidMpc+x)
        # apply the galaxy velocity correction
        z += (vel[0]/3e5)*(1+z)

#        haloId = str(int(haloId)) 
#        # find the plotting parameters based on the fof group   
#        if groups[haloId] == maxs[0]:
#            firstMassive = np.append(firstMassive, z)
#        elif groups[haloId] == maxs[1]:
#            secondMassive = np.append(secondMassive, z)
#        else:
#            lessMassive = np.append(lessMassive, z)

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
# - input -
#  d - distance to the object in Mpc
# - output -
#  z - the redshift of the object
def find_redshift(d, OmegaM=0.308, OmegaL=0.692):
    # this will help measure how long this step of the code will take
    starttime = time.time()
    OmegaK = 1-OmegaM-OmegaL
    # step to increment the redshift by
    dz = 0.001
    # initial redshift
    z = 3.05 
    # initialize the calculated distance to zero
    dCalculated = 0
    # while the calculated distance is less than the actual distance:
    while (dCalculated<d): 
        # increment the z that were testing for
        z += dz
        # create the array of zs for integration
        zarr = np.linspace(0,z, 200)
       
        # this is the integrand to get the comoving distance at some redshift 
        Ez = 1./np.sqrt(OmegaM*(1+zarr)**3+OmegaK*(1+zarr)**2+OmegaL)

        # this is the final distance
        # DH() is the hubble distance 
        dCalculated = cmToMpc*DH()*scipy.integrate.simps(Ez, zarr)

    # when the first step in z reached the required distance, return that value of z
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

    zarr = np.linspace(3.05, 3.12, 200)
#    plt.plot(zarr, twoGauss(zarr, *param1), 'k', linewidth=2)
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




# this is the field of influence for finding the column density of galaxies
def gal_kernel(x, y, x0, y0):

    rsqr = (x-x0)**2+(y-y0)**2
    return np.exp(-rsqr/10.)




# this will take in the positions that are rotated to the viewing angle
#  and return the list of indices that are with a fov of LRIS that is 
#  centered on the highest density region
def get_fov(positions):
    
    # size of the field of view
    width = 10.45
    height = 14.44

    # create a plain scatter plot of all of the halos
    plt.subplot(223)
    plt.xlim([-30,30])
    plt.ylim([-30,30])
    plt.xlabel(r"$H^{-1} \rm \ Mpc$", fontsize=16)
    plt.ylabel(r"$H^{-1} \rm \ Mpc$", fontsize=16)

  
    for pos in positions:
        plt.plot(pos[1], pos[2], 'ko', markersize=3)



    # create a density map of all of the galaxies
    Nx = 20
    Ny = 20
    x = np.linspace(-30,30, Nx)
    y = np.linspace(-30,30, Ny)
    X, Y = np.meshgrid(x, y)

    density = np.zeros((Nx, Ny), dtype='float')
    
    for ii in range(Nx):
        for jj in range(Ny):
            for pos in positions:
                density[ii,jj] += gal_kernel(x[ii], y[jj], pos[1], pos[2])

    plt.subplot(224)
    density = density[::-1,::-1]
    plt.contourf(X, Y, density, 100, cmap=plt.get_cmap('plasma'))
    plt.colorbar().set_label("Arbitrary Density", fontsize=16)
    plt.xlabel(r"$H^{-1} \rm \ Mpc$", fontsize=16)
    plt.ylabel(r"$H^{-1} \rm \ Mpc$", fontsize=16)


    # find the location of the peak density
    maxDensity = np.unravel_index(density.argmax(), density.shape)    
    print("maxDensity: {}".format(maxDensity))
    maxXY = [x[maxDensity[1]], y[maxDensity[0]]]
    print("maxXY: {}".format(maxXY))

    # the positions that will be returned
    subPositions = np.array([])
    for pos in positions:
        # check if it is in the right x locaiton
        if (pos[1] < maxXY[0]+width/2.) and (pos[1] > maxXY[0]-width/2.):
            # check if it is in the right y location
            if (pos[2] < maxXY[1]+height/2.) and (pos[2] > maxXY[1]-width/2.):
                subPositions = np.append(subPositions, pos)


    # plot the location of the max on the density plot
    plt.plot(maxXY[0], maxXY[1], 'w^', markersize=6)


    # reshape the subPositions array
    subPositions = np.reshape(subPositions, (-1,3))

    plt.subplot(223)
    for pos in subPositions:
        plt.plot(pos[1], pos[2], 'ro', markersize=4)

    
    np.set_printoptions(precision=2)
    return subPositions





if __name__=="__main__":
    # this will split up each of the individual descendant halos, and willl result
    #  in the random numbers being re-generated
#    split_query("./halos/MillenniumSQL_Full.dat")
#    find_fof("./halos/MillenniumSQL_Full.dat")    
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
