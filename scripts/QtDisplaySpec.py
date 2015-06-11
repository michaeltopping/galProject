import sys
import os
import matplotlib
import random
matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
from PyQt5.QtWidgets import *
import pyfits
import matplotlib.pyplot as plt
import time

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from math import *
from scipy.optimize import curve_fit
import smoothSpec
from scipy import interpolate

#constants
LyA = 1215.67

#logfile
logfile = open("QtDisplaySpec.log", "w")

#define a gaussian to fit to
def gauss(wavelengths, a, center, sigma):
	return a*np.exp(-(wavelengths-center)**2/(2*sigma**2))


#dataset class
class dataset():
	def __init__(self, folder, listFile, **kwargs):
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
		self.emRedshifts = np.array([])
		self.absRedshifts = np.array([])
		self.apertures = np.array([])
		for ii in range(np.shape(data)[0]):
			
			if not('usecols' in kwargs):
				data = np.genfromtxt(specListFile, dtype=np.str_)
				
# 				self.objects = np.append(self.objects, data[ii][2:-1])
# 				self.filenames = np.append(self.filenames, data[ii][2:-1])
				self.objects = np.append(self.objects, data[ii])
				self.filenames = np.append(self.filenames, data[ii])
				self.emRedshifts = np.append(self.emRedshifts, -2)
				self.absRedshifts = np.append(self.absRedshifts, -2)
				self.apertures = np.append(self.apertures, 1)
			else:
				data = np.genfromtxt(specListFile, dtype=np.str_, usecols=kwargs['usecols'], filling_values="xxx")

# 				fields = str(data[ii][0][2:-1])
# 				objects = str(data[ii][1][2:-1])
				fields = str(data[ii][0])
				objects = str(data[ii][1])
				try:
# 					emRedshift = float(data[ii][2][2:-1])
# 					absRedshift = float(data[ii][3][2:-1])
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
					
				
				self.objects = np.append(self.objects, objects)
				self.filenames = np.append(self.filenames, str(fields+".b."+objects+".msdc_v.fits"))
				self.emRedshifts = np.append(self.emRedshifts, emRedshift)
				self.absRedshifts = np.append(self.absRedshifts, absRedshift)
				self.apertures = np.append(self.apertures, aperture)

	
	
		
		self.spectrum = spectrum(self.folder, self.filenames[0], self.emRedshifts[0], self.absRedshifts[0], self.apertures[0])
			
		
	def setSpectrum(self, index):
		self.spectrum = spectrum(self.folder, self.filenames[index], self.emRedshifts[index], self.absRedshifts[0], self.apertures[index])
		


#spectrum class
class spectrum():
	def __init__(self,loc,  filename, emRedshift, absRedshift, aperture):
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
			lambdaMin = 3100.
			dlambda = 1.43
			Nlambda = np.shape(self.spec)[0]
		#Create an array of wavelengths from fits header wavelengths solution
		print("Creating wavelength solution")
		self.wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)

		
		#set limits on the wavelength that we are looking at
		minWavelength = 4500
		maxWavelength = 5300
		minWavelengthIndex = np.where(abs(self.wavelengths-minWavelength) == min(abs(self.wavelengths-minWavelength)))[0]
		maxWavelengthIndex = np.where(abs(self.wavelengths-maxWavelength) == min(abs(self.wavelengths-maxWavelength)))[0]

		print("Only using wavelegnth data from indices: ", minWavelengthIndex, " to ", maxWavelengthIndex)
		self.wavelengths = self.wavelengths[minWavelengthIndex:maxWavelengthIndex]
		self.spec = self.spec[minWavelengthIndex:maxWavelengthIndex]
		print(np.shape(self.spec), np.shape(self.wavelengths))

		
		self.getSmoothedSpectrum(51)
		#self.removeSpec(5600)
		#self.findRms(3)
		#self.removeDupPeaks()

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
		
		
	def fitLine(self, peak, plotBool):
		#how many wavelength points away from center we want to consider
		subSpecN = 10
		#create subarrays from center +/- subSpecN points
		
		peak = peak[0]
		print("Peak in fit: ", peak)
		subSpectra = self.spec[peak-subSpecN:peak+subSpecN]
		subWavelengths = self.wavelengths[peak-subSpecN:peak+subSpecN]

		#try to fit to gaussian
		#first, check to see if the line is close to the edge
		if (peak > subSpecN):
			maxSpec = max(subSpectra)

			fitParam, fitCov = curve_fit(gauss, subWavelengths-subWavelengths[subSpecN-1], subSpectra/maxSpec)
			#print out fit data
			print("A = ", fitParam[0])
			print("Center = ", fitParam[1]+subWavelengths[subSpecN-1], "Angstroms")
			print("Sigma = ", fitParam[2])
			z = (fitParam[1]+subWavelengths[subSpecN-1])/LyA - 1
			
			#create a plot of each spectrum and the fitted line.
			if (plotBool):
				print("Checking to see if directory exists")
				if not(os.path.isdir("./images/"+self.filename[0:-5])):
					os.makedirs("./images/"+self.filename[0:-5])
					print("Created directory")
				
				plt.cla()
				plt.plot(self.wavelengths, self.spec)
				plt.savefig("./images/"+self.filename[0:-5]+"/fullspec.png")
				plt.cla()
				plt.plot(subWavelengths, subSpectra/maxSpec)
				plt.plot(subWavelengths, fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[subSpecN-1]))**2/(2*fitParam[2]**2)))
				plt.savefig("./images/"+self.filename[0:-5]+"/linespec.png")
				
				if not(os.path.isdir("./images/all/")):
					os.makedirs("./images/all/")
				plt.cla()
				plt.plot(self.wavelengths, self.spec)
				plt.savefig("./images/all/"+str(z)+"_"+self.filename[0:-5]+"_fullspec.png")
				plt.cla()
				plt.plot(subWavelengths, subSpectra/maxSpec)
				plt.plot(subWavelengths, fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[subSpecN-1]))**2/(2*fitParam[2]**2)))
				plt.savefig("./images/all/"+str(z)+"_"+self.filename[0:-5]+"_linespec.png")
				


			return (fitParam[1]+subWavelengths[subSpecN-1])/LyA-1

			
		else:
			return -2



#define a base canvas class
class PlotCanvas(FigureCanvas):
	def __init__(self, parent=None, width=5, height=4, dpi = 100):
		fig = Figure(figsize=(width, height), dpi = dpi)
		self.axes = fig.add_subplot(111)
		
		
		FigureCanvas.__init__(self, fig)
		self.setParent(parent)
		
		FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)
	
	#for inheritance purposes
	def compute_initial_figure(self):
		pass
		
		

		

# #define a dynamic canvas class
class DynamicPlotCanvas(PlotCanvas):
	def __init__(self, spectrum,  *args, **kwargs):
		#get ready to make the plot
		#self.filename = filename
		#print(self.filename)
		PlotCanvas.__init__(self, *args, **kwargs)
		
		#create a clock that will refresh the plot every frame
		timer = QtCore.QTimer(self)
		timer.timeout.connect(self.update_figure)
		timer.start(50)
		
		#initialize variables
		self.spectrum = spectrum #create the spectrum object


		
	#update the figure, this is called efery frame
	def update_figure(self):

		#clear the plot
		self.axes.cla()
		
		#draw the spectrum
		self.axes.plot(self.spectrum.wavelengths, self.spectrum.spec, linewidth=1)
		minY = min(self.spectrum.spec)
		maxY = max(self.spectrum.spec)
		

		#if there are peaks listed in the peaks list, draw them
		if self.spectrum.peaks.size:
			for peak in self.spectrum.peaks:
				self.axes.plot([self.spectrum.wavelengths[peak], self.spectrum.wavelengths[peak]], [minY-10, maxY+10], 'k--', linewidth = 2)
		#draw rough values for positions of lines
		if self.spectrum.emPeaks.size:
			self.axes.plot([self.spectrum.emPeaks, self.spectrum.emPeaks], [0, maxY+10], 'r--', linewidth = 2)
		if self.spectrum.absPeaks.size:
			for peak in self.spectrum.absPeaks:
				self.axes.plot([peak, peak], [minY-10, 0], 'r--', linewidth = 2)
		
		#draw threshold spectrum
# 		if self.spectrum.threshSpec.size:
# 			self.axes.plot(self.spectrum.wavelengths, self.spectrum.threshSpec, 'k')

		#plotting window parameters
		self.axes.set_xlim([min(self.spectrum.wavelengths), max(self.spectrum.wavelengths)])
		self.axes.set_xlabel("Wavelengths [Angstroms]")
		self.axes.set_ylabel("Counts")
		self.axes.set_ylim([minY, maxY])
		#draw the figure
		self.draw()
		

		
		
		
#build the application window class
class ApplicationWindow(QMainWindow):
	def __init__(self):
		QMainWindow.__init__(self)
		#application attributes
		self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
		self.setWindowTitle("Line Finder")
		self.dataSetNum = 0
		
		#read in list of spectra from all.cat file
		specListFile = "../spec/notes/all.cat"
		self.specNum = 0
		data = np.loadtxt(specListFile, dtype=np.str_, usecols=(0,1, 2, 3))
		self.filenames = np.array([])
		self.emRedshifts = np.array([])
		self.absRedshifts = np.array([])
		for ii in range(np.shape(data)[0]):
			fields = str(data[ii][0][2:-1])
			objects = str(data[ii][1][2:-1])
			try:
				emRedshift = float(data[ii][2][2:-1])
				absRedshift = float(data[ii][3][2:-1])
			except:
				print("Error reading in redshift data")
				emRedshift = -2.0
				absRedshift = -2.0
				
			self.filenames = np.append(self.filenames, str(fields+".b."+objects+".msdc_v.fits"))
			self.emRedshifts = np.append(self.emRedshifts, emRedshift)
			self.absRedshifts = np.append(self.absRedshifts, absRedshift)
	#def __init__(self, folder, listFile, usecols, identifier=None):

		self.dataSet = [dataset("../spec/", "all.cat", usecols =(0, 1, 2, 3, 4)), 
						dataset("../shapley2003_spec/", "all.cat"),
						dataset("../shapley2006_spec/", "all.cat")]
		
		#more application attributes
		self.file_menu = QMenu('&File', self)
		self.file_menu.addAction('&Quit', self.fileQuit, QtCore.Qt.META + QtCore.Qt.Key_W)
		self.menuBar().addMenu(self.file_menu)

		
		#create the three buttons
		self.findSpecButton = QPushButton("Find Emission Lines")
		self.getRedshiftButton = QPushButton("Find z")
		self.nextSpecButton = QPushButton("Next Spectrum")
		self.prevSpecButton = QPushButton("Prev Spectrum")	
		self.createPlotsCheck = QCheckBox("Plots")
		self.getAllZButton = QPushButton("Find All Redshifts")
		self.loadDataButton = QPushButton("Load Spectra From File")
		self.cycleDataSetButton = QPushButton("Cycle DataSet")
		#create the main widget	
		self.main_widget = QWidget(self)
		
		#add the elements to the layout
		self.layout = QVBoxLayout(self.main_widget)

		self.specSelectLayout = QHBoxLayout()
		self.specSelectLayout.addWidget(self.prevSpecButton)
		self.specSelectLayout.addStretch(1)
		self.specSelectLayout.addWidget(self.findSpecButton)
		self.specSelectLayout.addStretch(1)
		self.specSelectLayout.addWidget(self.nextSpecButton)

		self.layout.addLayout(self.specSelectLayout)
		
		self.findZLayout = QHBoxLayout()
		self.Ztext = QLabel("<--- Finding Redshift Tools --->")
		self.Ztext.setAlignment(QtCore.Qt.AlignCenter)
		self.findZLayout.addWidget(self.getRedshiftButton)
		self.findZLayout.addStretch(1)

		self.findZLayout.addWidget(self.createPlotsCheck)
		self.findZLayout.addWidget(self.getAllZButton)

		self.layout.addLayout(self.findZLayout)
		
		
		self.layout.addWidget(self.loadDataButton)
		self.layout.addWidget(self.cycleDataSetButton)
		
		self.dynamic = DynamicPlotCanvas(self.dataSet[self.dataSetNum].spectrum, self.main_widget, width=10, height=4, dpi = 100)
		
		#self.layout.addWidget(self.getRedshiftButton)
		#self.layout.addWidget(self.getAllZButton)

		self.layout.addWidget(self.dynamic) #this is the spectrum window
		
		self.layout.addStretch(1)

		

		#connect the buttons so they can be clicked
		self.findSpecButton.clicked.connect(self.findLines)
		self.getRedshiftButton.clicked.connect(self.findZ)
		self.nextSpecButton.clicked.connect(self.nextSpec)
		self.prevSpecButton.clicked.connect(self.prevSpec)
		self.getAllZButton.clicked.connect(self.findAllZ)
		self.loadDataButton.clicked.connect(self.loadData)
		self.cycleDataSetButton.clicked.connect(self.cycleDataSet)
		
		#make the 
		self.main_widget.setFocus()
		self.setCentralWidget(self.main_widget)
		
	#quit the application
	def fileQuit(self):
		self.close()
		
	#also quit the application
	def closeEvent(self, ce):
		self.fileQuit()

	#find all the lines and plot them in the plotting window
	#findRms takes parameters of rms threshold for finding lines
	def findLines(self):
		self.dynamic.spectrum.findRms(4)
		self.dynamic.spectrum.removeDupPeaks()
		
	#find the redshift of current spectrum
	def findZ(self):
		#if no peaks have been searched for yet
		if (not self.dynamic.spectrum.peaks.size):
			self.dynamic.spectrum.findRms(4)
			self.dynamic.spectrum.removeDupPeaks()

		#if there are any peaks to be found
		if self.dynamic.spectrum.peaks.size:
			#get some local variables for data
			wavelengths = self.dynamic.spectrum.wavelengths
			spec = self.dynamic.spectrum.spec
			peaks = np.array(self.dynamic.spectrum.peaks)
			#if there are more than one peak, only look at the highest one
			highestPeak = np.where(spec == np.max(spec[peaks.astype(int)]))
			#fit the peak to a gaussian, return the center
			peakfit = self.dynamic.spectrum.fitLine(highestPeak)
			print(peakfit)

			
	#find redshifts of all spectra
	def findAllZ(self):
		starttime = time.time()
		lines = np.array([1260.42, 1303.3, 1334.5, 1393.76, 1402.77, 1526.7, 1549.5])
		W = np.array([1.63, 2.2, 1.72, 1.83, 0.81, 1.72, 3.03])
		f = np.array([1.007, 0.04887, .1278, .5140, .2553, .130, .1])
		plotBool = self.createPlotsCheck.isChecked()
		#empty array that will hold the redshifts
		zs = np.array([])
		laezs = np.array([])
		abszs = np.array([])
		objects = np.array([])
		nabszs = 0
		for dataset in self.dataSet:
		
			print(dataset)
			#loop through each file
			for file in range(len(dataset.filenames)):
				#create a local spectrum object for each file
				spec = spectrum(dataset.folder, dataset.filenames[file], dataset.emRedshifts[file], dataset.absRedshifts[file], dataset.apertures[file])
				#find all peaks in the file
				spec.findRms(4)
				spec.removeDupPeaks()
			
				#define local data
				wavelengths = spec.wavelengths
				data = spec.spec
				peaks = np.array(spec.peaks)
				#make sure we dont have any repeats.
				if not(dataset.objects[file] in objects):
					objects = np.append(objects, dataset.objects[file])
					#if there are any peaks in the spectrum
					if peaks.size:
						try:
							#pick the highest peak
							highestPeak = np.where(data == np.max(data[peaks.astype(int)]))
							#fit the peak to a gaussian and return the center
							peakfit = spec.fitLine(highestPeak, plotBool)
							print(peakfit)
							zs = np.append(zs, peakfit)
							#check if the objects is labeled as a Lyman alpha emitter
							if ('lae' in dataset.filenames[file]):
								laezs = np.append(laezs, peakfit)
							logfile.write(str(file)+ "	"+ dataset.filenames[file]+"	z=" +str(peakfit)+ '\n')
						except:
							print("Error finding peak")
					#if no peaks were found, use estimation from emission lines
					elif dataset.emRedshifts[file] > 0:
						try:
							#same fitting technique as before
							emRedshiftIndex = np.argmin(np.abs(wavelengths -(dataset.emRedshifts[file]+1)*LyA))
							print("Emredshift:",emRedshiftIndex, dataset.emRedshifts[file])
						
							peakfit = spec.fitLine([emRedshiftIndex], plotBool)
							zs = np.append(zs, peakfit)
							if ('lae' in dataset.filenames[file]):
								laezs = np.append(laezs, peakfit)
							logfile.write(str(file)+ "	"+ dataset.filenames[file]+"	z=" +str(peakfit)+ '\n')

						except RuntimeError:
							print("Unable to fit Gaussian")
							logfile.write(str(file) +"	"+ dataset.filenames[file]+ "	Unable to fit Gaussian"+ '\n')
					elif dataset.absRedshifts[file] > 0:
# 						try:
						print("Calculating absorption redshift")
						nabszs += 1
						zguess = dataset.absRedshifts[file]
						zlist = np.linspace(zguess*0.998,zguess*1.002, 1000)
						integrals = np.array([])
						for z in zlist:
							lineSpec = np.zeros(np.shape(spec.spec)[0])
							for kk in range(lines.size):
								line = ((z+1)*lines[kk])
								gaussian = f[kk]*np.exp(-((wavelengths - line)**2)/(2*W[kk]**2))
								lineSpec += gaussian
							multSpec = (spec.spec)*lineSpec
							integrals = np.append(integrals, np.trapz(wavelengths, multSpec))
						abszs = np.append(abszs, zlist[integrals.argmax()])
						logfile.write(str(file)+ "	"+ dataset.filenames[file]+" abs	z=" +str(zlist[integrals.argmax()])+ '\n')
						
						# except:
# 							print("Unable to fit Gaussian")
# 							logfile.write(str(file) +"	"+ dataset.filenames[file]+ "	Unable to calculate absorption redshift"+ '\n')

				else:
					print("This object is a duplicate")
					logfile.write(str(file)+ "	"+ dataset.filenames[file]+" Duplicate "+ '\n')

		print("Got ", zs.size, " Redshifts: ",zs)
		print("Got ", laezs.size, "LAE Redshifts: ",laezs)
		print("Got ", nabszs, "Absorption Redshifts: ", abszs)

		print("Found ", zs.size, " redshifts in ", time.time()-starttime, " seconds")
		
		#plot histogram
		#plt.hist(zs, bins=500, range=(2, 3.5),histtype="stepfilled",color='black' )

		plt.hist([zs, abszs], bins=500, range=(2, 3.5),histtype="stepfilled",color=['midnightblue', 'blue'],  stacked=True )
		#plt.ylim([0,20])
		#plt.hist(laezs, bins=500, range=(2, 3.5))
		#plt.hist(abszs, bins=500, range=(2, 3.5))
		plt.xlabel("z")
		plt.ylabel("Number")
		plt.xlim([3.0,3.2])
		plt.show()
			

	#switch data sets
	def cycleDataSet(self):
		self.dataSetNum += 1
		self.specNum = 0
		self.layout.removeWidget(self.dynamic)
		self.dynamic.deleteLater()
		self.dynamic = None
		self.dataSet[self.dataSetNum%3].setSpectrum(self.specNum)
		self.dynamic = DynamicPlotCanvas(self.dataSet[self.dataSetNum%3].spectrum, self.main_widget, width=10, height=4, dpi = 100)
		self.layout.addWidget(self.dynamic)

			
		
		
	#display the previous spectrum in the list
	def prevSpec(self):
		if self.specNum > 0:
			self.specNum -= 1
			#remove the current object
			self.layout.removeWidget(self.dynamic)
			#delete widget, free up memory
			self.dynamic.deleteLater()
			self.dynamic = None
			#create new spectrum pass
			self.dataSet[self.dataSetNum%3].setSpectrum(self.specNum)
			self.dynamic = DynamicPlotCanvas(self.dataSet[self.dataSetNum%3].spectrum, self.main_widget, width=10, height=4, dpi = 100)
			self.layout.addWidget(self.dynamic)
		else:
			pass
		
	#display the next spectrum in the list
	def nextSpec(self):
		self.specNum += 1
		#same as "prevSpec()" ^
		self.layout.removeWidget(self.dynamic)
		self.dynamic.deleteLater()
		self.dynamic = None
		self.dataSet[self.dataSetNum%3].setSpectrum(self.specNum)
		self.dynamic = DynamicPlotCanvas(self.dataSet[self.dataSetNum%3].spectrum, self.main_widget, width=10, height=4, dpi = 100)
		self.layout.addWidget(self.dynamic)
		
	def loadData(self):
		fname = QFileDialog.getOpenFileName( self, 'Select Files', '', "", "", QFileDialog.DontUseNativeDialog )
		print(fname[0])
		spec, header = pyfits.getdata(fname[0], 0, header=True)
		print(spec)

#main
if __name__ == '__main__':
	#start the application
	app = QApplication(sys.argv)
	
	aw = ApplicationWindow()
	aw.show()
	
	app.exec_()
	logfile.close()