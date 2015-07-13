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
from datasetClass import dataset
from spectrumClass import spectrum




#constants
LyA = 1215.67

#logfile
logfile = open("QtDisplaySpec.log", "w")

#define a gaussian to fit to
def gauss(wavelengths, a, center, sigma):
	return a*np.exp(-(wavelengths-center)**2/(2*sigma**2))


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

		#slider default
		sliderDefault = 30

		#read in list of spectra from all.cat file
		specListFile = "../spec/notes/all.cat"
		self.specNum = 0
		data = np.loadtxt(specListFile, dtype=np.str_, usecols=(0,1, 2, 3))
		self.filenames = np.array([])
		self.guessemRedshifts = np.array([])
		self.guessabsRedshifts = np.array([])
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
			self.guessemRedshifts = np.append(self.guessemRedshifts, emRedshift)
			self.guessabsRedshifts = np.append(self.guessabsRedshifts, absRedshift)
	#def __init__(self, folder, listFile, usecols, identifier=None):

		self.dataSet = [dataset("../spec/", "all.cat", "msdfc_v.fits", usecols =(0, 1, 2, 3, 4)),
						dataset("../shapley2003_spec/", "all.cat", "msdfcv.fits", usecols =(0, 1, 2, 3, 4)),
						dataset("../shapley2006_spec/", "all.cat", "br.fits", usecols =(0, 1, 2, 3, 4))]



		#more application attributes
		self.file_menu = QMenu('&File', self)
		self.file_menu.addAction('&Quit', self.fileQuit, QtCore.Qt.META + QtCore.Qt.Key_W)
		self.menuBar().addMenu(self.file_menu)


		#create the buttons
		self.findSpecButton = QPushButton("Find Emission Lines")
		self.getRedshiftButton = QPushButton("Find z")
		self.nextSpecButton = QPushButton("Next Spectrum")
		self.prevSpecButton = QPushButton("Prev Spectrum")
		self.createEmPlotsCheck = QCheckBox("Em Plots")
		self.createAbsPlotsCheck = QCheckBox("Abs Plots")
		self.getAllZButton = QPushButton("Find All Redshifts")
		self.loadDataButton = QPushButton("Load Spectra From File")
		self.cycleDataSetButton = QPushButton("Cycle DataSet")
		self.coordsButton = QPushButton("Display Coordinates")
		self.subSpecNSlider = QSlider(QtCore.Qt.Horizontal, self)
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
		self.sliderLabel=QLabel("SubSpecN: "+str(sliderDefault))
		self.subSpecNSlider.setMinimum(1)
		self.subSpecNSlider.setMaximum(50)

		self.subSpecNSlider.setSliderPosition(sliderDefault)
		#self.Ztext = QLabel("<--- Finding Redshift Tools --->")
		#self.Ztext.setAlignment(QtCore.Qt.AlignCenter)
		self.findZLayout.addWidget(self.getRedshiftButton)
		self.findZLayout.addStretch(1)
				
		self.findZLayout.addWidget(self.sliderLabel)
		self.findZLayout.addWidget(self.subSpecNSlider)
		self.findZLayout.addWidget(self.createEmPlotsCheck)
		self.findZLayout.addWidget(self.createAbsPlotsCheck)
		self.findZLayout.addWidget(self.getAllZButton)


		self.layout.addLayout(self.findZLayout)


		self.layout.addWidget(self.loadDataButton)
		self.layout.addWidget(self.cycleDataSetButton)
		self.layout.addWidget(self.coordsButton)

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
		self.subSpecNSlider.valueChanged[int].connect(self.changeSliderValue)
		self.coordsButton.clicked.connect(self.displayCoords)

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

	#activates when the value of the slider changes and gets the new value as input.
	#changes the label to the new value
	def changeSliderValue(self, value):
		self.sliderLabel.setText("SubSpecN: "+str(value))


	#find redshifts of all spectra
	def findAllZ(self):
		starttime = time.time()
		lines = np.array([1260.42, 1303.3, 1334.5, 1393.76, 1402.77, 1526.7, 1549.5])
		W = np.array([1.63, 2.2, 1.72, 1.83, 0.81, 1.72, 3.03])
		f = np.array([1.007, 0.04887, .1278, .5140, .2553, .130, .1])
		emPlotBool = self.createEmPlotsCheck.isChecked()
		absPlotBool = self.createAbsPlotsCheck.isChecked()
		#empty array that will hold the redshifts
		zs = np.array([])
		laezs = np.array([])
		abszs = np.array([])
		objects = np.array([])
		nabszs = 0
		subSpecN = self.subSpecNSlider.value()
		for dataset in self.dataSet:

			print(dataset)
			#loop through each file
			for file in range(len(dataset.filenames)):
				#create a local spectrum object for each file
				spec = spectrum(dataset.folder, dataset.filenames[file], dataset.guessemRedshifts[file], dataset.guessabsRedshifts[file], dataset.apertures[file])
				#find all peaks in the file
				spec.findRms(4)
				spec.removeDupPeaks()

				#define local data
				wavelengths = spec.wavelengths
				data = spec.spec
				peaks = np.array(spec.peaks)
				#make sure we dont have any repeats.
				#if not(dataset.objects[file] in objects):
				objects = np.append(objects, dataset.objects[file])
				#if there are any peaks in the spectrum
				if peaks.size:
					try:
						#pick the highest peak
						highestPeak = np.where(data == np.max(data[peaks.astype(int)]))
						#fit the peak to a gaussian and return the center
						peakfit = spec.fitLine(highestPeak, emPlotBool, subSpecN)
						print(peakfit)


						#check if the line is double peaked, then it will have different behavior when calculating systematic redshift
						doublePeaks = spec.isDoublePeaked(spec.wavelengths, spec.spec, spec.wavelengths[highestPeak])
						if len(doublePeaks) == 2:
							doublePeaked = True
						else:
							doublePeaked = False

						if doublePeaked:
							trough = spec.spec[min(doublePeaks):max(doublePeaks)].argmin()+min(doublePeaks)
							peakfit = spec.wavelengths[trough]/LyA - 1



						zs = np.append(zs, peakfit)
						#check if the objects is labeled as a Lyman alpha emitter
						if ('lae' in dataset.filenames[file]):
							laezs = np.append(laezs, peakfit)
							dataset.galaxies[dataset.objects[file]].addType("lae")
						logfile.write(str(file)+ "	"+ dataset.filenames[file]+"	z=" +str(peakfit)+ '\n')
						#absz = spec.fitAbsLine(highestPeak, plotBool)
						#abszs = np.append(abszs, absz)
						dataset.galaxies[dataset.objects[file]].addRedshift("em", peakfit)
						#dataset.galaxies[dataset.objects[file]].addRedshift("abs", absz)
					except:
						print("Error finding peak")
				#if no peaks were found, use estimation from emission lines
				elif dataset.guessemRedshifts[file] > 0:
					try:
						#same fitting technique as before
						emRedshiftIndex = np.argmin(np.abs(wavelengths -(dataset.guessemRedshifts[file]+1)*LyA))
						print("Emredshift:",emRedshiftIndex, dataset.guessemRedshifts[file])

						peakfit = spec.fitLine([emRedshiftIndex], emPlotBool, subSpecN)

						zs = np.append(zs, peakfit)
						if ('lae' in dataset.filenames[file]):
							laezs = np.append(laezs, peakfit)
						logfile.write(str(file)+ "	"+ dataset.filenames[file]+"	z=" +str(peakfit)+ '\n')

						dataset.galaxies[dataset.objects[file]].addRedshift("em", peakfit)

					except RuntimeError:
						print("Unable to fit Gaussian")
						logfile.write(str(file) +"	"+ dataset.filenames[file]+ "	Unable to fit Gaussian"+ '\n')
				if dataset.guessabsRedshifts[file] > 0:
# 						try:
					print("Calculating absorption redshift")
					nabszs += 1
					zguess = dataset.guessabsRedshifts[file]

					absz = spec.fitAbsLine(zguess, absPlotBool)

					abszs = np.append(abszs, absz)
					logfile.write(str(file)+ "	"+ dataset.filenames[file]+" abs	z=" +str(absz)+ '\n')
					dataset.galaxies[dataset.objects[file]].addRedshift("abs", absz)



		print("Got ", zs.size, " Redshifts: ",zs)
		print("Got ", laezs.size, "LAE Redshifts: ",laezs)
		print("Got ", abszs.size, "Absorption Redshifts: ", abszs)

		print("Found ", zs.size, " redshifts in ", time.time()-starttime, " seconds")

		self.plotHistogram()

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

	def printData(self):
		for dataset in self.dataSet:
			dataset.outputData()

	def plotHistogram(self):
		zs = np.array([])
		abszs = np.array([])

		for dataset in self.dataSet:
			for gal in dataset.galaxies:
				dataset.galaxies[gal].systematicShift()
				z = np.average(dataset.galaxies[gal].sysRedshift)
				absz = np.average(dataset.galaxies[gal].absRedshifts)
				if isnan(z):
					z = 0
				if isnan(absz):
					absz = 0
				if z>0:
					print("object:", gal, "of type", dataset.galaxies[gal].type, "has redshifts:", dataset.galaxies[gal].emRedshifts, "with average:", z)
				zs = np.append(zs, z)
				abszs = np.append(abszs, absz)
		#plot histogram
		plt.cla()
		#for emission and absorption
		#plt.hist([zs, abszs], bins=500, range=(2, 3.5),histtype="step",color=['blue', 'red'],  stacked=False, linewidth=3 )
		#for emission only
		plt.hist(zs, bins=500, range=(2, 3.5), histtype="stepfilled" ,color='blue')
		plt.xlabel("z")
		plt.ylabel("Number")
		plt.xlim([3.03,3.12])
		plt.show()


	def displayCoords(self):
		#add positional information to datasets
		self.dataSet[0].setAllCoords("../spec/mask_design/")

#main
if __name__ == '__main__':
	#start the application
	app = QApplication(sys.argv)

	aw = ApplicationWindow()
	aw.show()

	app.exec_()
	logfile.close()
