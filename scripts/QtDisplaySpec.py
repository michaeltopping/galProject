import sys
import matplotlib
import random
matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
from PyQt5.QtWidgets import *
import pyfits

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from math import *
from scipy.optimize import curve_fit
import smoothSpec
from scipy import interpolate

#spectrum class
class spectrum():
	def __init__(self, filename, emRedshift, absRedshift):
		self.emRedshift = emRedshift
		self.absRedshift = absRedshift
		print(emRedshift)
		#read in the spectrum data
		#create an empty list of redshifts
		zlist = np.array([])
		self.peaks = np.array([])
		self.emPeaks = np.array([])
		self.absPeaks = np.array([])
		self.threshSpec = np.array([])
		
		#Read in data into spec, and header data into header
		print("Reading in data from "+filename)
		self.spec, header = pyfits.getdata("../spec/"+filename, 0, header=True)
		if (float(header['NAXIS']) > 1):
			nSpectra = np.shape(self.spec)
		else:
			nSpectra = 1
			
		if nSpectra == 1:
			self.spec = self.spec
		else:
			self.spec = self.spec[0]
		print("Nspectra: ", nSpectra)
		print("Read in", nSpectra, "spectra.")
		
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
		minWavelength = 4000
		maxWavelength = 5400
		minWavelengthIndex = np.where(abs(self.wavelengths-minWavelength) == min(abs(self.wavelengths-minWavelength)))[0]
		maxWavelengthIndex = np.where(abs(self.wavelengths-maxWavelength) == min(abs(self.wavelengths-maxWavelength)))[0]

		print("Only using wavelegnth data from indices: ", minWavelengthIndex, " to ", maxWavelengthIndex)
		self.wavelengths = self.wavelengths[minWavelengthIndex:maxWavelengthIndex]
		self.spec = self.spec[minWavelengthIndex:maxWavelengthIndex]
		
		self.getSmoothedSpectrum(101)
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

		self.emPeaks = (self.emRedshift+1)*1216.
		self.absPeaks = (self.absRedshift+1)*np.array([1303.3,1334.5])

		
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
		print("Adding peak to list")
		if jj != 0:
			peakIndices = np.append(peakIndices, int((self.peaks[jj-1]+self.peaks[jj-1-nInPeak])/2))
		else:
			print("No peaks to add to list")
		print("Peak added to list")
		
		self.peaks = (peakIndices)
		print("Peaks: ", self.peaks)


#define a base canvas class
class PlotCanvas(FigureCanvas):
	def __init__(self, parent=None, width=5, height=4, dpi = 100):
		fig = Figure(figsize=(width, height), dpi = dpi)
		self.axes = fig.add_subplot(111)
		#redraw axes when plot is called
		#self.axes.hold(False)
		
		
		FigureCanvas.__init__(self, fig)
		self.setParent(parent)
		
		FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
		FigureCanvas.updateGeometry(self)
	
	def compute_initial_figure(self):
		pass
		
		
#define a static canvas class
class StaticPlotCanvas(PlotCanvas):
	def compute_initial_figure(self):
		x = np.arange(0, 100, 0.01)
		y = np.sin(x)


		self.axes.plot(x, y)
		

# #define a dynamic canvas class
class DynamicPlotCanvas(PlotCanvas):
	def __init__(self, filename, emZ, absZ,  *args, **kwargs):
		#get ready to make the plot
		self.filename = filename
		print(self.filename)
		PlotCanvas.__init__(self, *args, **kwargs)
		
		#create a clock that will refresh the plot every frame
		timer = QtCore.QTimer(self)
		timer.timeout.connect(self.update_figure)
		timer.start(50)
		
		#initialize variables
		self.time = 0
		self.freq = 0
		self.spectrum = spectrum(self.filename, emZ, absZ) #create the spectrum object


		
	#update the figure, this is called efery frame
	def update_figure(self):
		#increase a time variable
		self.time += 1
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
		if self.spectrum.threshSpec.size:
			print(np.shape(self.spectrum.wavelengths), np.shape(self.spectrum.threshSpec))
			self.axes.plot(self.spectrum.wavelengths, self.spectrum.threshSpec, 'k')

		#plotting window parameters
		self.axes.set_xlim([min(self.spectrum.wavelengths), max(self.spectrum.wavelengths)])
		self.axes.set_xlabel("Wavelengths [Angstroms]")
		self.axes.set_ylabel("Counts")
		self.axes.set_ylim([minY-10, maxY+10])
		#draw the figure
		self.draw()
		

		
		
		
#build the application window class
class ApplicationWindow(QMainWindow):
	def __init__(self):
		QMainWindow.__init__(self)
		#application attributes
		self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
		self.setWindowTitle("Line Finder")
		
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
			
		print(self.emRedshifts, self.absRedshifts)
		
		#more application attributes
		self.file_menu = QMenu('&File', self)
		self.file_menu.addAction('&Quit', self.fileQuit, QtCore.Qt.META + QtCore.Qt.Key_W)
		self.menuBar().addMenu(self.file_menu)
		
		#create the three buttons
		self.findSpecButton = QPushButton("Find Emission Lines")
		self.nextSpecButton = QPushButton("Next Spectrum")
		self.prevSpecButton = QPushButton("Previous Spectrum")	
		#create the main widget	
		self.main_widget = QWidget(self)
		
		#add the elements to the layout
		self.layout = QVBoxLayout(self.main_widget)
		self.dynamic = DynamicPlotCanvas(str(self.filenames[0]), self.emRedshifts[0], self.absRedshifts[0], self.main_widget, width=10, height=4, dpi = 100)
		self.layout.addWidget(self.findSpecButton)
		self.layout.addWidget(self.nextSpecButton)
		self.layout.addWidget(self.prevSpecButton)
		self.layout.addWidget(self.dynamic) #this is the spectrum window

		#connect the buttons so they can be clicked
		self.findSpecButton.clicked.connect(self.findLines)
		self.nextSpecButton.clicked.connect(self.nextSpec)
		self.prevSpecButton.clicked.connect(self.prevSpec)
		
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
		self.dynamic.spectrum.findRms(3)
		self.dynamic.spectrum.removeDupPeaks()
		
	#display the previous spectrum in the list
	def prevSpec(self):
		if self.specNum > 0:
			self.specNum -= 1
			self.layout.removeWidget(self.dynamic)
			self.dynamic.deleteLater()
			self.dynamic = None
			self.dynamic = DynamicPlotCanvas(str(self.filenames[self.specNum]), self.emRedshifts[self.specNum], self.absRedshifts[self.specNum], self.main_widget, width=10, height=4, dpi = 100)
			self.layout.addWidget(self.dynamic)
		else:
			pass
		
	#display the next spectrum in the list
	def nextSpec(self):
		self.specNum += 1
		self.layout.removeWidget(self.dynamic)
		self.dynamic.deleteLater()
		self.dynamic = None
		self.dynamic = DynamicPlotCanvas(str(self.filenames[self.specNum]), self.emRedshifts[self.specNum], self.absRedshifts[self.specNum], self.main_widget, width=10, height=4, dpi = 100)
		self.layout.addWidget(self.dynamic)

#main
if __name__ == '__main__':
	#start the application
	app = QApplication(sys.argv)
	
	aw = ApplicationWindow()
	aw.show()
	
	app.exec_()