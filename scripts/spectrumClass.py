import numpy as np
import matplotlib.pyplot as plt
import pyfits
from math import *
from scipy.optimize import curve_fit
import smoothSpec
from scipy import interpolate
import os
import sys

#constants
LyA = 1215.67

#absorption lines
lines = np.array([1260.42, 1303.3, 1334.5, 1393.76, 1402.77, 1526.7, 1549.5])
W = np.array([1.63, 2.2, 1.72, 1.83, 0.81, 1.72, 3.03])
f = np.array([1.007, 0.04887, .1278, .5140, .2553, .130, .1])

#define a gaussian to fit to
def gauss(wavelengths, a, center, sigma, b, c):
	return c*wavelengths+b+a*np.exp(-(wavelengths-center)**2/(2*sigma**2))


#spectrum class
#procedures:
# init-
#  parameters:
#   loc - directory that the spectrum is in
#   filename - filename of spectrum file
#   emRedshift - guess for the redshift based on emission features
#   absRedshift - guess of the redshift based on absorption features
#   
class spectrum():
	def __init__(self, loc,  filename, emRedshift, absRedshift, aperture):
		self.emRedshift = emRedshift
		self.absRedshift = absRedshift
		print(emRedshift)
		self.filename = filename
		#read in the spectrum data
		#create an empty list of redshifts
		zlist = np.array([])
		self.peaks = np.array([])
		self.emPeaks = np.array([])
		self.absPeaks = np.array([])
		self.threshSpec = np.array([])
		
		#Read in data into spec, and header data into header
		print("Reading in data from "+filename)
		self.spec, header = pyfits.getdata(loc+filename, 0, header=True)
		print(np.shape(self.spec))
		if (float(header['NAXIS']) > 1):
			nSpectra = np.shape(self.spec)[0]
		else:
			nSpectra = 1
		print(nSpectra, " Spectra detected")
		if nSpectra == 1:
			self.spec = self.spec
		else:
			self.spec = self.spec[aperture-1]
		print("Nspectra: ", nSpectra)
		print("Read in", nSpectra, "spectra.")
		print("Shape of spectra: ", np.shape(self.spec))
		
		self.spec[np.where(self.spec > 1000)] = 0
		self.spec[np.where(self.spec < -1000)] = 0
	

		
		
		#Read in and split up wavelength solution information in the 'WAT2_002' fits header
		try:
			lambdaMin = float(header['CRVAL1'])
			dlambda = float(header['CD1_1'])
			Nlambda = float(header['NAXIS1'])
		except:
			#alternate wavelength solution uses WAT2_002 header parameter
			foundMin = False
			index = 0
			while not(foundMin):
				index += 1
				try:
					print("One element: ", float(header['WAT2_002'].split()[index]))
					if (float(header['WAT2_002'].split()[index]) > 2000.0):
						foundMin = True
						print("Min Wavelength is: ", float(header['WAT2_002'].split()[index]))
					
				except:
					pass
					
		
			lambdaMin = float(header['WAT2_002'].split()[index])
			dlambda = float(header['WAT2_002'].split()[index+1])
			Nlambda = np.shape(self.spec)[0]
		#Create an array of wavelengths from fits header wavelengths solution
		print("Creating wavelength solution")
		self.wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)


		#we need the full spectrum and wavelength data for some things
		self.fullspec = self.spec
		self.fullwavelengths = self.wavelengths

		
		#set limits on the wavelength that we are looking at
		minWavelength = 4500
		maxWavelength = 5300
		minWavelengthIndex = np.where(abs(self.wavelengths-minWavelength) == min(abs(self.wavelengths-minWavelength)))[0]
		maxWavelengthIndex = np.where(abs(self.wavelengths-maxWavelength) == min(abs(self.wavelengths-maxWavelength)))[0]

		print("Only using wavelegnth data from indices: ", minWavelengthIndex, " to ", maxWavelengthIndex)
		self.wavelengths = self.wavelengths[minWavelengthIndex:maxWavelengthIndex]
		self.spec = self.spec[minWavelengthIndex:maxWavelengthIndex]

		
		self.getSmoothedSpectrum(51)

	#LyAguess is in angstroms
	def isDoublePeaked(self, wavelengths, spec, LyAguess):
		dspec = np.array([])
		dlambda = np.array([])

		#set limits on the wavelength that we are looking at
		minWavelength = LyAguess-10
		maxWavelength = LyAguess+10
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
	

		


		print("Npeaks: ", len(list(np.where(peakSpec[peaks]>0.8)[0])))
		print("Calculated peaks at: ", peakWavelengths[peaks])
		LyApeaks = []
		for peak in list(np.where(peakSpec[peaks]>0.8)[0]):
			peak = peaks[peak]
			print("Trying ", peakWavelengths[peak])
			print("Found peaks at: ", wavelengths[np.where(wavelengths == peakWavelengths[peak])[0][0]])
			LyApeaks.append(np.where(wavelengths == peakWavelengths[peak])[0][0])
		print("Done finding peaks --------------")
		return LyApeaks
# 		if len(list(np.where(spec[peaks]>0.5)[0])) == 2:
# 			return True
# 		else:
# 			return False


	def getSmoothedSpectrum(self, smoothFac):
		#get a smoothed spectrum
		
		interp = interpolate.splrep(smoothSpec.smooth(self.wavelengths, smoothFac), smoothSpec.smooth(self.spec, smoothFac))
		lambdaNew = np.linspace(min(self.wavelengths), max(self.wavelengths), len(self.spec))
		self.smoothSpec = interpolate.splev(lambdaNew, interp, der=0)


	def removeSpec(self, center):
		print("Removing atmosphere line")
		#find index of dichroic boundary
		dichroIndex = np.where(abs(self.wavelengths-center)==min(abs(self.wavelengths-center)))[0]
		#remove this boundary from the spectrum
		maskSize = 30
		for index in np.arange(dichroIndex-maskSize,dichroIndex+maskSize, 1):
			self.spec[index] = (self.spec[dichroIndex]+self.spec[dichroIndex+2*maskSize])/2.	


	def findRms(self, rmsThresh):
		#find the rms deviation in the spectra
		rmsarr = np.array([])
		self.threshSpec = np.array([])
		rms = np.std(self.spec - self.smoothSpec)
		print("RMS error in spectra", 0, "is", rms)
		#pick out data points that are 5x above the rms value
		#find the rms as a function of wavelength
		Nslice = 5
		self.peaks = np.array([])
		sliceSize = int(np.shape(self.spec)[0]/Nslice)
		print("Size of slice: ", sliceSize)
		for ii in range(Nslice):
			rms = np.std(self.spec[ii*sliceSize:(ii+1)*sliceSize] - self.smoothSpec[ii*sliceSize:(ii+1)*sliceSize])
			rmsarr = np.append(rmsarr, rms)
			#print("Adding peaks ", (ii*sliceSize)+np.where(self.spec[0][ii*sliceSize:(ii+1)*sliceSize] - self.smoothSpec[ii*sliceSize:(ii+1)*sliceSize] > rmsThresh*rms)[0])
			self.peaks = np.append(self.peaks, (ii*sliceSize)+np.where(self.spec[ii*sliceSize:(ii+1)*sliceSize] - self.smoothSpec[ii*sliceSize:(ii+1)*sliceSize] > rmsThresh*rms)[0])

			self.threshSpec = np.append(self.threshSpec, self.smoothSpec[ii*sliceSize:(ii+1)*sliceSize+1]+rmsThresh*rmsarr[ii])

		self.threshSpec = self.threshSpec[0:-1]

		self.emPeaks = (self.emRedshift+1)*LyA
		self.absPeaks = (self.absRedshift+1)*np.array([1303.3,1334.5])
		print("Found peaks at ", self.peaks)

		
	def removeDupPeaks(self):
		lastpeak = 0
		nInPeak = 0
		peakIndices = np.array([])
		#find index near center of lines picked out of lise of 'peaks' array
		jj = 0
		for jj in range(np.size(self.peaks)):
			print("Looking at peak number ", jj, "-------------------------")
			#if this is the first value we are trying
			if lastpeak == 0:
				lastpeak = self.wavelengths[self.peaks[jj]]
				print("This is the first peak")
			#if this point is in the same spectral line as the last point
			elif ( (self.wavelengths[self.peaks[jj]]-lastpeak) < 10):
				lastpeak = self.wavelengths[self.peaks[jj]]
				nInPeak += 1
				print("This is close to the last peak")
				print("Current peak: ", self.wavelengths[self.peaks[jj]])
				print("Last peak: ", lastpeak)
			#if this point is in a new spectral line
			elif ( (self.wavelengths[self.peaks[jj]]-lastpeak) >= 10):
				print("Added new peak, distance since last peak is: ", (self.wavelengths[self.peaks[jj]]-lastpeak))
				print("Now on peak jj: ", self.peaks[jj], " last peak was ", self.peaks[jj-1], " first index in the peak is ", self.peaks[jj-nInPeak-1])
				lastpeak = self.wavelengths[self.peaks[jj]]
				#add the previous point to a list of line centers
				peakIndices = np.append(peakIndices, int((self.peaks[jj-1]+self.peaks[jj-1-nInPeak])/2))
				print("Number of indices in the peak: ", nInPeak)
				print("Now on peak number ", jj)
				nInPeak = 0
				print("Added peak: ", int((self.peaks[jj-1]+self.peaks[jj-1-nInPeak])/2))
	
		#required to add the last line to the list of centers
		print("Adding last peak to list after", nInPeak, " peaks")
		if jj != 0:
			peakIndices = np.append(peakIndices, int((self.peaks[jj]+self.peaks[jj-nInPeak])/2))
		else:
			print("No peaks to add to list")
		#print("Peak added to list", int((self.peaks[jj]+self.peaks[jj-nInPeak])/2))
		
		self.peaks = (peakIndices)
		print("Reduced peaks at peaks: ", self.peaks)
		
		
	def fitLine(self, peak, plotBool, subSpecN):
		#how many wavelength points away from center we want to consider
		#subSpecN = 15
		#create subarrays from center +/- subSpecN points
		
		peak = peak[0]
		print("Peak in fit: ", peak)
		subSpectra = self.spec[peak-subSpecN:peak+subSpecN]
		subWavelengths = self.wavelengths[peak-subSpecN:peak+subSpecN]
		plotSpectra = self.spec[peak-3*subSpecN:peak+3*subSpecN]
		plotWavelengths = self.wavelengths[peak-3*subSpecN:peak+3*subSpecN]





		doublePeaks = self.isDoublePeaked(self.wavelengths, self.spec, self.wavelengths[peak])
		if len(doublePeaks) == 2:
			doublePeaked = True
		else:
			doublePeaked = False
			
		print("Peaks in LyA: ", doublePeaks)
		print("Corresponds to: ", self.wavelengths[doublePeaks])

		if doublePeaked:
			trough = self.spec[min(doublePeaks):max(doublePeaks)].argmin()+min(doublePeaks)
			

		#try to fit to gaussian
		#first, check to see if the line is close to the edge
		if (peak > subSpecN):
			maxSpec = max(subSpectra)

			fitParam, fitCov = curve_fit(gauss, (subWavelengths-subWavelengths[subSpecN-1]), subSpectra/maxSpec)
			#print out fit data
			print("A = ", fitParam[0])
			print("Center = ", fitParam[1]+subWavelengths[subSpecN-1], "Angstroms")
			print("Sigma = ", fitParam[2])
			print("Continuum = ", fitParam[3])
			print("Continuum slope = ", fitParam[4])
			z = (fitParam[1]+subWavelengths[subSpecN-1])/LyA - 1
			
			#create a plot of each spectrum and the fitted line.
			if (plotBool):
				print("Checking to see if directory exists")
				if not(os.path.isdir("./images/"+self.filename[0:-5])):
					os.makedirs("./images/"+self.filename[0:-5])
					print("Created directory")
				if fitParam[0]>0:
					plt.cla()
					plt.plot(self.wavelengths, self.spec)
					plt.savefig("./images/"+self.filename[0:-5]+"/fullspec.png")
					plt.cla()
					plt.plot(subWavelengths, subSpectra/maxSpec, linewidth=2)
					plt.plot(subWavelengths, (subWavelengths-subWavelengths[subSpecN-1])*fitParam[4]+fitParam[3]+fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[subSpecN-1]))**2/(2*fitParam[2]**2)), linewidth=2, label="z=%.3f" % z+"d"*doublePeaked)
					plt.legend()
					
					plt.savefig("./images/"+self.filename[0:-5]+"/linespec.png")
				
					if not(os.path.isdir("./images/all/")):
						os.makedirs("./images/all/")
					plt.cla()
					plt.plot(self.wavelengths, self.spec)
					plt.savefig("./images/all/"+str(z)+"_"+self.filename[0:-5]+"_fullspec.png")
					plt.cla()
					plt.plot(subWavelengths, subSpectra/maxSpec, linewidth=2)
					plt.plot(subWavelengths, (subWavelengths-subWavelengths[subSpecN-1])*fitParam[4]+fitParam[3]+fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[subSpecN-1]))**2/(2*fitParam[2]**2)), linewidth=2, label="z=%.3f" % z+"d"*doublePeaked)
					plt.plot(subWavelengths, (subWavelengths-subWavelengths[subSpecN-1])*fitParam[4]+fitParam[3]+fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[subSpecN-1]))**2/(2*fitParam[2]**2)) - subSpectra/maxSpec, 'r')
					if doublePeaked:
						plt.plot([self.wavelengths[trough], self.wavelengths[trough]], [0,1], 'k')
					plt.legend()
					if doublePeaked:
						plt.plot(self.wavelengths[doublePeaks], self.spec[doublePeaks]/max(subSpectra), 'ko')
					plt.savefig("./images/all/"+str(z)+"d"*doublePeaked+"_"+self.filename[0:-5]+"_linespec.png")
				
			if not doublePeaked:
				if fitParam[0]<0:
					return -1*(fitParam[1]+subWavelengths[subSpecN-1])/LyA+1
				elif fitParam[0]>0:
					return (fitParam[1]+subWavelengths[subSpecN-1])/LyA-1
			else:
				return self.wavelengths[trough]/LyA - 1

			
		else:
			return -2
	
	def fitAbsLine(self, peak, plotBool):
		print("Calculating absorption redshift")
		zguess = self.wavelengths[peak[0]]/LyA-1
		print("Zguess is ", zguess)
		zlist = np.linspace(zguess*0.998,zguess*1.002, 100)
		integrals = np.array([])
		zstep = 0
		for z in zlist:
			zstep += 1
			sys.stdout.write("\rCalculating absorption redshift: |"+"-"*int(zstep/4)+"."*(25-int(zstep/4))+"|%d%%" % int(zstep))
			sys.stdout.flush()
			lineSpec = np.zeros(np.shape(self.fullspec)[0])
			for kk in range(lines.size):
				line = ((z+1)*lines[kk])
				gaussian = f[kk]*np.exp(-((self.fullwavelengths - line)**2)/(2*W[kk]**2))
				lineSpec += gaussian
			multSpec = (self.fullspec)*lineSpec
			integrals = np.append(integrals, np.trapz(self.fullwavelengths, multSpec))
			
		absZ = zlist[integrals.argmax()]
			
			
		if (plotBool):
			print("Checking to see if directory exists")
			if not(os.path.isdir("./images/"+self.filename[0:-5])):
				os.makedirs("./images/"+self.filename[0:-5])
				print("Created directory")
			
			plt.cla()
			plt.plot(self.fullwavelengths, self.fullspec)
			for peak in lines*(absZ+1):
				plt.plot([peak, peak], [min(self.fullspec)-10, 0], 'r--', linewidth = 2)
			plt.savefig("./images/"+self.filename[0:-5]+"/absorbspec.png")
			
			if not(os.path.isdir("./images/all/absorption/")):
				os.makedirs("./images/all/absorption/")
			plt.cla()
			plt.plot(self.fullwavelengths, self.fullspec)
			for peak in lines*(absZ+1):
				plt.plot([peak, peak], [min(self.fullspec)-10, 0], 'r--', linewidth = 2)
			plt.savefig("./images/all/absorption/"+str(absZ)+"_"+self.filename[0:-5]+"_absorbspec.png")

			
			
			
		return absZ






