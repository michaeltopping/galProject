import numpy as np
import matplotlib.pyplot as plt








#define class to hold galaxy data
class galaxy():
    def __init__(self, name, RA, dec, z, emZ, absZ, LyAtype):
        self.name = name
        self.RA = RA
        self.dec = dec
        self.z = z
        self.emZ = emZ
        self.absZ = absZ
        self.LyAtype = LyAtype

def dataInit():
    galaxies={}
    #read in positional information
    posData = np.genfromtxt("./masterPositions.dat", dtype=str)

    #read in redshift data
    zData = np.genfromtxt("table.dat", dtype=str, skip_header=1)

    #loop through each galaxy that we have redshift information for and create a new gala    xy object
    for ii in range(np.shape(zData)[0]):
        #the name of the galaxy
        zName = zData[ii][0]
        #this is the systemic redshift
        z = float(zData[ii][1])
        #this is em and abs redshift data
        emZ = float(zData[ii][2])
        absZ = float(zData[ii][3])
        #if there is no redshift information, z will be -9.999
        if z > 0:
            #the index in the position data array that contains this galaxy
            posIndex = np.where(zName == posData[:,0])[0]
            #check if the galaxy exists in the position array
            if posIndex:
                #create a new galaxy object in the galaxies dictionary
                galaxies[zName] = galaxy(zName, float(posData[posIndex,1][0]), float(posData[posIndex,2][0]), z, emZ, absZ, zData[ii][4])
    print("Done adding galaxy data to {} galaxies.".format(len(galaxies)))
    return galaxies



def plot_all_hist():
    # create a dictionary of all the galaxy data
    galaxies = dataInit()

    # create an array that contains all of the galaxy redshifts
    redshifts = np.array([])
    # loop through the galaxy dictionary
    for gal in galaxies:
        if galaxies[gal].z > 0:
            redshifts = np.append(redshifts, galaxies[gal].z)


    # do the plotting of the redshift here
    plt.hist(redshifts, bins=30, range=(3.05, 3.12), histtype="stepfilled", color="blue")
    plt.xlabel("z")
    plt.ylabel("N")
    plt.show()


def LAE_LBG_hist():
    # create a dictionary of all the galaxy data
    galaxies = dataInit()

    # create an array for the LBG and LAE redshifts
    LAEzs = np.array([])
    LBGzs = np.array([])


    # loop through galaxies in the dictionary and separate the redshifts
    #  into LAE and LBG redshifts
    for gal in galaxies:
        # first check if an LBG by the galaxy ID
        if gal.startswith("C") or gal.startswith("D") or gal.startswith("M"):
            # this object is an LBG
            LBGzs = np.append(LBGzs, galaxies[gal].z)
        elif gal.startswith("NB") or gal.startswith("lae") or gal[0].isdigit():
            # this object is an LAE
            LAEzs = np.append(LAEzs, galaxies[gal].z)
        else:
            print("Error, object {} is not LAE or LBG".format(gal))


    # plot the histogram here
    plt.hist([LBGzs, LAEzs], bins=30, range=(3.05, 3.12), stacked=True, histtype="stepfilled", label=["LBGs", "LAEs"])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("N")
    plt.show()
     
    


def LyA_shape_hist():
   # create a dictionary of all galaxy data
    galaxies = dataInit() 
    
    #create lists for single and double peaked LyA
    singlePeakZs = np.array([])
    doublePeakZs = np.array([])
    
    # loop through all of the galaxies
    for gal in galaxies:
        #check if the galaxy is double peaked
        if galaxies[gal].LyAtype == "d":
            # galaxy had double peaked Lyman Alpha
            doublePeakZs = np.append(doublePeakZs, galaxies[gal].z) 
        else:
            # galaxy does not have double peaked Lyman Alpha
            singlePeakZs = np.append(singlePeakZs, galaxies[gal].z)


    # plot the histogram here
    plt.hist([doublePeakZs, singlePeakZs], bins=30, range=(3.05, 3.12), stacked=True, histtype="stepfilled", label=["Double Peaked LyA", "Single Peaked LyA"])
    plt.legend(loc="upper left")
    plt.xlabel("z")
    plt.ylabel("N")
    plt.show()




def abs_em_hist():
    # create a dictionary of all galaxy data
    galaxies = dataInit()

    # create arrays for galaxies with emission only, absorption only, and both
    absZs = np.array([])
    emZs = np.array([])
    bothZs = np.array([])
    
    # loop through galaxies
    for gal in galaxies:
        #check if the galaxy has both em and abs redshift
        if ( galaxies[gal].emZ > 0 ) and ( galaxies[gal].absZ > 0 ):
            # this galaxy has both em and abs redshifts
            bothZs = np.append(bothZs, galaxies[gal].z)
        # check if the galaxy only has em Redshift
        elif galaxies[gal].emZ > 0:
            # this galaxy only has an emission redshfits    
            emZs = np.append(emZs, galaxies[gal].z)
        # check if the galaxy only has abs redshift 
        elif galaxies[gal].absZ > 0:
            # this galaxy only has absorption redshifts 
            absZs = np.append(absZs, galaxies[gal].z)
        # if the galaxy has no redshifts
        else:
            print("Galaxy {} has no redshift measurements.".format(gal)) 


    # plot the histogram here
    plt.hist([absZs, bothZs, emZs], bins=30, range=(3.05, 3.12), stacked=True, histtype="stepfilled", label=["z_abs only", "z_em and z_abs", "z_em only"])
    plt.legend(loc="upper left")
    plt.xlabel("z")
    plt.ylabel("N")
    plt.show()

if __name__ == "__main__":
   abs_em_hist()
else:
    print("Plotting options:")
    print("  plot_all_hist()")
    print("  LAE_LBG_hist()")
    print("  LyA_shape_hist()")
    print("  abs_em_hist()")
 


