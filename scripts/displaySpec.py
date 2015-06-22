#----------------------------------
#
#This script will read in a fits file, then plot it
#
#----------------------------------
import numpy as np
import pylab as plt
import pyfits
from smoothSpec import smooth
from scipy import interpolate
import random

def isDoublePeaked(wavelengths, spec, LyAguess):
	dspec = np.array([])
	dlambda = np.array([])

	#set limits on the wavelength that we are looking at
	minWavelength = LyAguess-30
	maxWavelength = LyAguess+30
	minWavelengthIndex = np.where(abs(wavelengths-minWavelength) == min(abs(wavelengths-minWavelength)))[0][0]
	maxWavelengthIndex = np.where(abs(wavelengths-maxWavelength) == min(abs(wavelengths-maxWavelength)))[0][0]


	peakWavelengths = wavelengths[minWavelengthIndex:maxWavelengthIndex]
	peakSpec = spec[minWavelengthIndex:maxWavelengthIndex]	
	peakSpec = peakSpec/max(peakSpec)
	for ii in range(len(peakWavelengths)-1):
		dspec = np.append(dspec, peakSpec[ii+1]-peakSpec[ii])
		dlambda = np.append(dlambda, (peakWavelengths[ii+1]+peakWavelengths[ii])/2.)
	negs=[]
	for jj in range(len(dspec)):
		if dspec[jj]<0:
			negs.append(jj)

	pos = [x+1 for x in list(set(range(len(spec)-1))-set(negs))]
	peaks = list(set(pos) & set(negs))


	print peaks


	print "Npeaks: ", len(list(np.where(peakSpec[peaks]>0.5)[0]))
	for peak in list(np.where(peakSpec[peaks]>0.5)[0]):
		print wavelengths[np.where(wavelengths == peakWavelengths[peak])[0]][0], peakWavelengths[peak]
		

	return list(np.where(spec[peaks]>0.5)[0])
	


#define the location of the data file
DIR = "../shapley2006_spec/"
filename = "ssa22.D17.br.fits"

#Read in data into spec, and header data into header
print "Reading in data from "+DIR+filename
spec, header = pyfits.getdata(DIR+filename, 0, header=True)
nSpectra = np.shape(spec)[0]
print "Read in", nSpectra, "spectra."

lambdaMin = float(header['CRVAL1'])
dlambda = float(header['CD1_1'])
Nlambda = float(header['NAXIS1'])

#Create an array of wavelengths from fits header wavelengths solution
print "Creating wavelength solution"
wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)
smoothFac = 101
interp = interpolate.splrep(smooth(wavelengths, smoothFac), smooth(spec, smoothFac))
lambdaNew = np.linspace(min(wavelengths), max(wavelengths), len(spec))
specNew = interpolate.splev(lambdaNew, interp, der=0)



#plot data
# for ii in range(nSpectra):
# 	print "Plotting spectrum number", ii
# 	plt.plot(wavelengths, spec[ii])
#plt.plot(wavelengths, spec)
# plt.plot(wavelengths, specNew)

#change plot limits
# plt.ylim([-50, 100])
# 
# #Show the plot
# plt.show()

lines = np.array([1260.42, 1303.3, 1334.5, 1393.76, 1402.77, 1526.7, 1549.5])
W = np.array([1.63, 2.2, 1.72, 1.83, 0.81, 1.72, 3.03])
f = np.array([1.007, 0.04887, .1278, .5140, .2553, .130, .1])
Nzs = 1000
zguess = 3.091
#zs = np.linspace(zguess*0.998,zguess*1.002, Nzs)
zs = np.linspace(zguess*0.75,zguess*1.25, Nzs)

lambdaMin = 3100.
dlambda = 1.43
integrals = np.array([])
totSpec = np.array([spec[0]])


waves = np.fft.fft(spec)
freq = np.fft.fftfreq(lambdaNew.shape[-1])
print freq[waves.argmin()]
# plt.plot(freq, waves.real)
# plt.show()

totspec = spec
totints = np.zeros(Nzs)
# N = 100
# for jj in range(N):
# 	std = np.std(spec)
# 	noiseSpec = np.array([])
# 	for ii in range(wavelengths.size):
# 		noiseSpec = np.append(noiseSpec, 2*std*(random.random()-0.5))
# 	
# 	totspec = totspec + spec+noiseSpec
# 	
# 	integrals = np.array([])
# 	for z in zs:
# 		lineSpec = np.zeros(np.shape(spec)[0])
# 		for kk in range(lines.size):
# 			line = ((z+1)*lines[kk])
# 			gaussian = f[kk]*np.exp(-((wavelengths - line)**2)/(2*(W[kk]/dlambda)**2))
# 			lineSpec += gaussian
# 		multSpec = (spec+noiseSpec-specNew)*lineSpec
# 		integrals = np.append(integrals, np.trapz(wavelengths, multSpec))
# 
# 
# 	totints = totints + integrals
# 	print "Finished calculating number ", jj
	
# plt.plot(wavelengths, totspec/(N+1))

bool = isDoublePeaked(wavelengths, spec, 4973.)
print bool

plt.plot(wavelengths, spec/max(spec))
plt.xlim([4500,5500])
plt.show()
# 
# integrals = np.array([])
# for z in zs:
# 	lineSpec = np.zeros(np.shape(spec)[0])
# 	for kk in range(lines.size):
# 		line = ((z+1)*lines[kk])
# 		gaussian = f[kk]*np.exp(-((wavelengths - line)**2)/(2*W[kk]**2))
# 		lineSpec += gaussian
# 	multSpec = (spec)*lineSpec
# 	integrals = np.append(integrals, np.trapz(wavelengths, multSpec))
# 
# #plt.plot(zs, totints/N)
# print "Redshift: ", zs[integrals.argmax()]
# plt.plot(zs, integrals)
# plt.show()

#plt.plot(zs, totints/N-integrals)
#plt.show()


# for ii in range(spec.size-1):
# 
# 	totSpec = np.append(totSpec, totSpec[ii]+spec[ii+1])
# 	
# 
# 	
# plt.plot(wavelengths, lineSpec)
# plt.show()
# plt.plot(wavelengths, totSpec)
# plt.show()
# plt.plot(zs, integrals)
# plt.show()
# print lineSpec



# length = wavelengths.size
# 
# der = np.array([])
# waveder = np.array([])
# for ii in range(length-2):
# 	der = np.append(der, spec[0][ii+1]-spec[0][ii])
# 	waveder = np.append(waveder, (wavelengths[ii+1]+wavelengths[ii])/2.)
# 	
# plt.plot(smooth(waveder, 3), smooth(der,3))
# plt.show()
	

