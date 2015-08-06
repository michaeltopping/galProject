import numpy as np
import matplotlib.pyplot as plt
from spectrumClass import spectrum
import pyfits
from math import *
from scipy.optimize import curve_fit
import smoothSpec
from scipy import interpolate
import os
from galaxy import Object
from matplotlib.colors import LogNorm

#constants
LyA = 1215.67

#define a gaussian to fit to
def gauss(wavelengths, a, center, sigma):
    return a*np.exp(-(wavelengths-center)**2/(2*sigma**2))

#dataset class
#this will be used if data from more than one place are being used, and we want to keep them separate
class dataset():
    def __init__(self, folder, listFile, suffix, **kwargs):
        self.folder = folder
        self.listFile = listFile
        specListFile = self.folder+self.listFile

        print("Reading in new dataset from", specListFile)

            
        #read in data from listFile, then make array of spectrum objects that big
        self.nSpectra = 10
        self.spectra = []

        #find the filename, emmission redshift, and absorption redshift from the listfile
    
        #read in list of spectra from all.cat file
        data = np.genfromtxt(specListFile, dtype=np.str_, usecols=(0,))
        self.filenames = np.array([])
        self.objects = np.array([])
        self.guessemRedshifts = np.array([])
        self.guessabsRedshifts = np.array([])
        self.emRedshifts = np.array([])
        self.absRedshifts = np.array([])
        self.apertures = np.array([])
        self.galaxies = {}
        objectlist = []

        #loop through each line in the "all.cat" file
        for ii in range(np.shape(data)[0]):
            
            #check the format of the file
            if not('usecols' in kwargs):
                data = np.genfromtxt(specListFile, dtype=np.str_)

                self.objects = np.append(self.objects, data[ii])
                self.filenames = np.append(self.filenames, data[ii])
                self.guessemRedshifts = np.append(self.guessemRedshifts, -2)
                self.guessabsRedshifts = np.append(self.guessabsRedshifts, -2)
                self.apertures = np.append(self.apertures, 1)

            #if it is in a 5 column format do this
            #the columns are:
            #mask - object - emission z - absorption z - notes
            else:
                data = np.genfromtxt(specListFile, dtype=np.str_, usecols=kwargs['usecols'], filling_values="xxx")

                object = str(data[ii][0])
                filename = str(data[ii][1])
                try:

                    emRedshift = float(data[ii][2])
                    absRedshift = float(data[ii][3])
                except:
                    print("Error reading in redshift data")
                    emRedshift = -2.0
                    absRedshift = -2.0
                try:
                    note = data[ii][4]
                except:
                    print("Error reading note")
                
                if "ap" in note:
                    try:    
                        aperture = int(note[2])
                    except IndexError:
                        aperture = 1
                else:
                    aperture = 1
                if not(object in objectlist):
                    objectlist.append(object)
                    self.galaxies[object] = Object()

                    
                #set member variables

                self.objects = np.append(self.objects, object)
                self.filenames = np.append(self.filenames, filename)
                self.guessemRedshifts = np.append(self.guessemRedshifts, emRedshift)
                self.guessabsRedshifts = np.append(self.guessabsRedshifts, absRedshift)
                self.apertures = np.append(self.apertures, aperture)
        #create a spectrum object as the first object in the list
        self.spectrum = spectrum(self.folder, self.filenames[0], self.objects[0], self.guessemRedshifts[0], self.guessabsRedshifts[0], self.apertures[0])
            
        
    #change the current spectrum to the spectrum identified with "index"
    def setSpectrum(self, index):
        self.spectrum = spectrum(self.folder, self.filenames[index], self.objects[index], self.guessemRedshifts[index], self.guessabsRedshifts[0], self.apertures[index])
        

    #this will print out the type of object each galaxy is, and its redshift.
    def outputData(self):
        for ii in range(self.filenames.size):
            print(self.filenames[ii]+" emission z="+str(self.galaxies[self.objects[ii]].emRedshifts) + " absorption z="+str(self.galaxies[self.objects[ii]].absRedshifts))



    #add positional information to the objects
    def setAllCoords(self, folder):

        #read in slit mask data
        for filename in os.listdir(folder):
            data = np.genfromtxt(folder+filename, dtype="str")
            for ii in range(np.shape(data)[0]):
                obj = data[ii][0]
                RA = float(data[ii][3])*15+float(data[ii][4])*15./60. + float(data[ii][5])*15/60./60.
                dec = float(data[ii][6])+float(data[ii][7])/60.+float(data[ii][8])/60./60.
                print("Coordinates for: {obj:11s} RA={RA:3f} dec={dec:3f}".format(obj=obj, RA=RA, dec=dec))
                if obj in self.galaxies:
                    self.galaxies[obj].setCoords(RA, dec)

        for gal in self.galaxies:
            x = (self.galaxies[gal].RA-334.3852)/(-1.388888e-5)+3100.
            y = (self.galaxies[gal].dec-0.19403)/(-1.388888e-5)+3100.
            plt.scatter(x, y, facecolors='none', edgecolors='r')
            #plt.text(self.galaxies[gal].RA, self.galaxies[gal].dec, gal)

        plt.xlabel("Right Ascension")
        plt.ylabel("Declination")
        # plt.xlim([334.32,334.43])
        # plt.ylim([0.175, 0.315])
        plt.xlim([0,6000])
        plt.ylim([0, 6000])

        #read in the hubble image
        image_data = pyfits.getdata("../image_data/HST_10405_07_ACS_WFC_F814W_drz.fits")
        image_data[image_data<0]=0

        image_data += 0.1

        plt.imshow(image_data, cmap='gray', vmax = 2*np.average(image_data), norm=LogNorm())

        plt.savefig("coords.eps", format="eps", dpi=1000)
        plt.show()
