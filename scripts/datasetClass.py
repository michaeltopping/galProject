import numpy as np
import matplotlib.pyplot as plt
from spectrumClass import spectrum
import pyfits
from math import *
from scipy.optimize import curve_fit
import smoothSpec
from scipy import interpolate
import os
from galaxy import galaxy

#constants
LyA = 1215.67

#define a gaussian to fit to
def gauss(wavelengths, a, center, sigma):
	return a*np.exp(-(wavelengths-center)**2/(2*sigma**2))

#dataset class
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
		for ii in range(np.shape(data)[0]):
			
			if not('usecols' in kwargs):
				data = np.genfromtxt(specListFile, dtype=np.str_)

				self.objects = np.append(self.objects, data[ii])
				self.filenames = np.append(self.filenames, data[ii])
				self.guessemRedshifts = np.append(self.guessemRedshifts, -2)
				self.guessabsRedshifts = np.append(self.guessabsRedshifts, -2)
				self.apertures = np.append(self.apertures, 1)
			else:
				data = np.genfromtxt(specListFile, dtype=np.str_, usecols=kwargs['usecols'], filling_values="xxx")

				fields = str(data[ii][0])
				object = str(data[ii][1])
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
					aperture = int(note[2])
				else:
					aperture = 1
				if not(object in objectlist):
					objectlist.append(object)
					self.galaxies[object] = galaxy()
					
				self.objects = np.append(self.objects, object)
				self.filenames = np.append(self.filenames, str(fields+"."+object+"."+suffix))
				self.guessemRedshifts = np.append(self.guessemRedshifts, emRedshift)
				self.guessabsRedshifts = np.append(self.guessabsRedshifts, absRedshift)
				self.apertures = np.append(self.apertures, aperture)
				
	
		self.spectrum = spectrum(self.folder, self.filenames[0], self.guessemRedshifts[0], self.guessabsRedshifts[0], self.apertures[0])
			
		
	def setSpectrum(self, index):
		self.spectrum = spectrum(self.folder, self.filenames[index], self.guessemRedshifts[index], self.guessabsRedshifts[0], self.apertures[index])
		

	def outputData(self):
		for ii in range(self.filenames.size):
			print(self.filenames[ii]+" emission z="+str(self.galaxies[self.objects[ii]].emRedshifts) + " absorption z="+str(self.galaxies[self.objects[ii]].absRedshifts))
