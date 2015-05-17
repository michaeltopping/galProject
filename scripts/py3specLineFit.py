#----------------------------------
#
#This script will read in a fits file, and fit all of the lines
#Also, plot the spectra
#----------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from math import *
from scipy.optimize import curve_fit
import smoothSpec
from scipy import interpolate



#define the location of the data file
DIR = "../spec/"
#filename = "ssaly_1.b.C12.msdc_v.fits"
# filename = "ssaly_1.b.79495.msdc_v.fits"
# filename = "test.fits"

#define a gaussian to fit to
def gauss(wavelengths, a, center, sigma):
	return a*np.exp(-(wavelengths-center)**2/(2*sigma**2))


#LyAlpha
#usage
#filename - name of the fits file that contains spectra
#rmsThresh - how many times the rms value of the spectra you would like to look
#plotBool - True/False would you like to plot the spectra
def LyAlpha(filename, rmsThresh, plotBool):

	#create an empty list of redshifts
	zlist = np.array([])

	#Read in data into spec, and header data into header
	print("Reading in data from "+DIR+filename)
	spec, header = pyfits.getdata(DIR+filename, 0, header=True)
	if (float(header['NAXIS']) > 1):
		nSpectra = np.shape(spec)[0]
	else:
		nSpectra = 1
	print("Nspectra: ", nSpectra)
	print("Read in", nSpectra, "spectra.")
	

	#Read in and split up wavelength solution information in the 'WAT2_002' fits header
	lambdaMin = float(header['CRVAL1'])
	dlambda = float(header['CD1_1'])
	Nlambda = float(header['NAXIS1'])

	#Create an array of wavelengths from fits header wavelengths solution
	print("Creating wavelength solution")
	wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)

	print(range(nSpectra))
	
	#plot/analyze data
	for ii in range(nSpectra):
		print(spec)
		if nSpectra == 1:
			spectrum = spec
		else:
			spectrum = spec[ii]
			
		print("Creating figure")
		#plt.figure(figsize=(13, 6))
		print("created figure")
		#change plot parameters
		#plt.xlabel("Wavelength $\AA$")
		#plt.ylabel("Counts")
		print("Plotting spectrum number", ii)


		sliceBool = True
		if sliceBool == True:
			#set limits on the wavelength that we are looking at
			minWavelength = 4000
			maxWavelength = 5400
			minWavelengthIndex = np.where(abs(wavelengths-minWavelength) == min(abs(wavelengths-minWavelength)))[0]
			maxWavelengthIndex = np.where(abs(wavelengths-maxWavelength) == min(abs(wavelengths-maxWavelength)))[0]

			print("Only using wavelegnth data from indices: ", minWavelengthIndex, " to ", maxWavelengthIndex)
			wavelengths = wavelengths[minWavelengthIndex:maxWavelengthIndex]
			spectrum = spectrum[minWavelengthIndex:maxWavelengthIndex]
	
		
		#get a smoothed spectrum
		smoothFac = 101
		interp = interpolate.splrep(smoothSpec.smooth(wavelengths, smoothFac), smoothSpec.smooth(spectrum, smoothFac))
		lambdaNew = np.linspace(min(wavelengths), max(wavelengths), len(spectrum))
		specNew = interpolate.splev(lambdaNew, interp, der=0)
		
		# #find index of dichroic boundary
# 		dichroIndex = np.where(abs(wavelengths-5600)==min(abs(wavelengths-5600)))[0]
# 		#remove this boundary from the spectrum
# 		maskSize = 30
# 		for index in np.arange(dichroIndex-maskSize,dichroIndex+maskSize, 1):
# 			spectrum[index] = (spectrum[dichroIndex]+spectrum[dichroIndex+2*maskSize])/2.	
# 	
		#find the rms deviation in the spectra
		rmsarr = np.array([])
		rms = np.std(spectrum - specNew)
		print("RMS error in spectra", ii, "is", rms)
		#pick out data points that are 5x above the rms value
		#find the rms as a function of wavelength
		Nslice = 3
		peaks = []
		sliceSize = int(np.shape(spectrum)[0]/Nslice)
		print("Size of slice: ", sliceSize)
		for ii in range(10):
			rms = np.std(spectrum[ii*sliceSize:(ii+1)*sliceSize] - specNew[ii*sliceSize:(ii+1)*sliceSize])
			rmsarr = np.append(rmsarr, rms)
			print("Adding peaks ", (ii*sliceSize)+np.where(spectrum[ii*sliceSize:(ii+1)*sliceSize] - specNew[ii*sliceSize:(ii+1)*sliceSize] > rmsThresh*rms)[0])
			peaks = np.append(peaks, (ii*sliceSize)+np.where(spectrum[ii*sliceSize:(ii+1)*sliceSize] - specNew[ii*sliceSize:(ii+1)*sliceSize] > rmsThresh*rms)[0])

		waverms = np.linspace(min(wavelengths), max(wavelengths), Nslice)

		print("There are ", len(peaks), " at indices ", dlambda*peaks+lambdaMin)
		lastpeak = 0
		nInPeak = 0
		peakIndices = np.array([])
		#find index near center of lines picked out of lise of 'peaks' array
		jj = 0
		for jj in range(np.size(peaks)):
			print("Looking at peak number ", jj, "-------------------------")
			#if this is the first value we are trying
			if lastpeak == 0:
				lastpeak = wavelengths[peaks[jj]]
				print("This is the first peak")
			#if this point is in the same spectral line as the last point
			elif ( (wavelengths[peaks[jj]]-lastpeak) < 10):
				lastpeak = wavelengths[peaks[jj]]
				nInPeak += 1
				print("This is close to the last peak")
				print("Current peak: ", wavelengths[peaks[jj]])
				print("Last peak: ", lastpeak)
			#if this point is in a new spectral line
			elif ( (wavelengths[peaks[jj]]-lastpeak) >= 10):
				print("Added new peak, distance since last peak is: ", (wavelengths[peaks[jj]]-lastpeak))
				print("Now on peak jj: ", peaks[jj], " last peak was ", peaks[jj-1], " first index in the peak is ", peaks[jj-nInPeak-1])
				lastpeak = wavelengths[peaks[jj]]
				#add the previous point to a list of line centers
				peakIndices = np.append(peakIndices, int((peaks[jj-1]+peaks[jj-1-nInPeak])/2))
				print("Number of indices in the peak: ", nInPeak)
				print("Now on peak number ", jj)
				nInPeak = 0
				print("Added peak: ", int((peaks[jj-1]+peaks[jj-1-nInPeak])/2))
	
		#required to add the last line to the list of centers
		print("Adding peak to list")
		if jj != 0:
			peakIndices = np.append(peakIndices, int((peaks[jj-1]+peaks[jj-1-nInPeak])/2))
		else:
			print("No peaks to add to list")
		print("Peak added to list")
	
		print("Reduced peak list: ", peakIndices*dlambda + lambdaMin)

		#plotting
		for peak in peakIndices:
			x = [wavelengths[peakIndices.astype(int)], wavelengths[peakIndices.astype(int)]]
			y = [-100, 500]
			#plt.plot(x,y, 'k--')

		#plt.plot(wavelengths, spectrum, 'k')
		#plt.plot(lambdaNew, specNew, 'r')

		#for rmsVal in range(Nslice):
			#plt.plot(lambdaNew[rmsVal*sliceSize:(rmsVal+1)*sliceSize], specNew[rmsVal*sliceSize:(rmsVal+1)*sliceSize]+rmsThresh*rmsarr[rmsVal], 'b')
			
			#print("Something")
			
		print("Finding max flux value")
		#find the max flux value
		#maxY = max(specNew[rmsVal*sliceSize:(rmsVal+1)*sliceSize]+rmsThresh*rmsarr[rmsVal])
		#add the filename to the plot
		#plt.text(minWavelength+10, maxY-.05*maxY, filename)
		#add which aperature is being plotted
		#plt.text(minWavelength+10, maxY-.1*maxY, "Aperature: "+str(ii+1))
		#plt.ylim([-50, maxY+10])
		#plt.xlim([minWavelength, maxWavelength])


		
	
		#fitting the lines to a gaussian
		npeak = 0
		
		for peak in peakIndices:	
			#current peak
			npeak = npeak +1
			
			#change plot parameters
			#plt.xlabel("Wavelength $\AA$")
			#plt.ylabel("Counts")
			
			#how many wavelength points away from center we want to consider
			subSpecN = 20
			#create subarrays from center +/- 20 points
			subSpectra = spectrum[peak-subSpecN:peak+subSpecN]
			subWavelengths = wavelengths[peak-subSpecN:peak+subSpecN]

			
			#try to fit to gaussian
			try:
				fitParam, fitCov = curve_fit(gauss, subWavelengths-subWavelengths[21], subSpectra)
				if ((fitParam[2] > 0) and (fitParam[2] < 10)):
					print("Creating new figure")
					#plt.figure()
					#plt.plot(subWavelengths, fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[21]))**2/(2*fitParam[2]**2)), label="Center: "+"{0:.2f}".format(fitParam[1]+subWavelengths[21])+"\n"+"z={0:.3f}".format((fitParam[1]+subWavelengths[21])/1216-1))
					#plt.legend(loc="upper left")
					#plot
					#plt.plot(subWavelengths, subSpectra, 'k')
					zlist = np.append(zlist, (fitParam[1]+subWavelengths[21])/1216-1)
					
				else:
					print("Was not able to create a good enough fit to a gaussian.")
				#print out fit data
				print("For peak number ", npeak, " out of ", len(peakIndices))
				print("A = ", fitParam[0])
				print("Center = ", fitParam[1]+subWavelengths[21], "Angstroms")
				print("Sigma = ", fitParam[2])
			except:
				print("Error fitting a gaussian")
		
		
		
		


	#change plot parameters
	#plt.xlabel("Wavelength $\AA$")
	#plt.ylabel("Counts")

	#Show the plot
	#if (plotBool == True):
	#	plt.show()
		
	return zlist
		
		
if __name__ == "__main__":
	filename = "../spec/notes/all.cat"
	data = np.loadtxt(filename, dtype=np.str_, usecols=(0,1))
	Zs = np.array([])
	plotBool = False
	nchecked = 0
	nfailed = 0
	#loop through all data
	for ii in range(np.shape(data)[0]):
		field = data[ii][0][2:-1]
		object = data[ii][1][2:-1]
		print(object, object[0:3])
		if (object[0:3] == "lae"):
			print("Going")
			filename = field+".b."+object+".msdc_v.fits"
			nchecked += 1
			try:
				print(filename+"----------------------------------------------")
				Zs = np.append(Zs, LyAlpha(filename, 4, plotBool))
			except:
				print("Something Went Wrong!")
				nfailed += 1
				
	print("List of redshifts based on LyA: ", Zs)
	print("Checked ", nchecked, " with ", nfailed, " failures")
	plt.hist(Zs, bins=40)
	plt.xlim([2.8, 3.3])
	plt.xlabel("z")
	plt.ylabel("N")
	plt.show()
