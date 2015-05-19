#----------------------------------
#
#This script will read in a fits file, then plot it
#
#----------------------------------
import numpy as np
import pylab as plt
import pyfits
from smoothSpec import smooth

#define the location of the data file
DIR = "../spec/"
filename = "ssaly_1.b.79495.msdc_v.fits"

#Read in data into spec, and header data into header
print "Reading in data from "+DIR+filename
spec, header = pyfits.getdata(DIR+filename, 0, header=True)
nSpectra = np.shape(spec)[0]
print "Read in", nSpectra, "spectra."

#Read in and split up wavelength solution information in the 'WAT2_002' fits header
waveSol = header['WAT2_002'].split()

#Create an array of wavelengths from fits header wavelengths solution
print "Creating wavelength solution"
wavelengths = np.arange(float(waveSol[4]), float(waveSol[4]) + float(waveSol[5])*float(waveSol[6]), float(waveSol[5]))



#plot data
for ii in range(nSpectra):
	print "Plotting spectrum number", ii
	plt.plot(wavelengths, spec[ii])


#change plot limits
plt.ylim([-50, 100])

#Show the plot
plt.show()

length = wavelengths.size

der = np.array([])
waveder = np.array([])
for ii in range(length-2):
	der = np.append(der, spec[0][ii+1]-spec[0][ii])
	waveder = np.append(waveder, (wavelengths[ii+1]+wavelengths[ii])/2.)
	
plt.plot(smooth(waveder, 3), smooth(der,3))
plt.show()
	

