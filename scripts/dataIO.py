import numpy as np








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
        self.R = 0

def galaxy_init():
    galaxies={}
    #read in positional information
    posData = np.genfromtxt("./masterPositions.dat", dtype=str)

    #read in redshift data
    zData = np.genfromtxt("table.dat", dtype=str, skip_header=1)

    #loop through each galaxy that we have redshift information for and create a new gala        xy object
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



# this will read in all of the photometry and return a dictionary with 
#  all of the galaxies that have photometry and their R-band Magnitudes
def read_photometry():
    # create the dictionary that will hold the magnitude data
    galaxyMagnitudes = {}
    # first read in all of the LBG magnitudes from the .current file
    filename = "../spec/mask_design/z.ssa22a.current"
    data = np.genfromtxt(filename, dtype=str, skip_header=1)
    # loop through all input data and assign magnitudes
    for row in data:
        if not row[2] in galaxyMagnitudes:
            galaxyMagnitudes[row[2]] = row[4]

   
    # now we will go through the objects in Nestor et al. 2011
    #  both LAEs and LBGs with and without NB detections
    #  will need to get some key to convert between their naming convention 
    filelist = ["../photometry/LAE_withNB3640.txt", "../photometry/LAE_withoutNB3640.txt", 
                "../photometry/LBG_withNB3640.txt", "../photometry/LBG_withoutNB3640.txt"]
    # loop through the files and add the galaxies to the dictionary
    for filename in filelist:
        print("Now looking at file: {}".format(filename))
        data = np.genfromtxt(filename, dtype=str, skip_header=5, usecols=(0, 6))
        
        # loop through all input data and assign magnitudes
        for row in data:
            if not row[0] in galaxyMagnitudes:
                galaxyMagnitudes[row[0]] = row[1] 
            # check if the object has already been read in
            elif row[0] in galaxyMagnitudes:
                print("Duplicate object found: {}".format(row[0]))

    # print out the dictionary
#    print(galaxyMagnitudes)






if __name__ =="__main__":
    read_photometry()

