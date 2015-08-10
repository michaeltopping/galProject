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
    def __init__(self, loc,  filename, objName, emRedshift, absRedshift, aperture):
        self.emRedshift = emRedshift
        self.absRedshift = absRedshift
        print(emRedshift)
        self.filename = filename
        self.objName = objName
        #read in the spectrum data
        #create an empty list of redshifts
        zlist = np.array([])
        self.peaks = np.array([])
        self.emPeaks = np.array([])
        self.absPeaks = np.array([])
        self.threshSpec = np.array([])
        

        #parameters for the fitting of LyA
        self.lineConst = 0
        self.lineLinear = 0
        self.lineCenter = 0
        self.lineSigma = 0
        self.lineAmplitude = 0

        #equivalent width
        self.EW = 0

        #type of peak
        self.peakType = ''

        #Read in data into spec, and header data into header
        print("Reading in data from "+filename)
        self.spec, header = pyfits.getdata(loc+filename, 0, header=True)

        #check how many apertures are in the image
        if (float(header['NAXIS']) > 1):
            nSpectra = np.shape(self.spec)[0]
        else:
            nSpectra = 1
        print(nSpectra, " Spectra detected")
        if nSpectra == 1:
            self.spec = self.spec
        else:
            self.spec = self.spec[aperture-1]

        #clip bad data
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
        minWavelengthIndex = np.where(abs(self.wavelengths-minWavelength) ==
                               min(abs(self.wavelengths-minWavelength)))[0]
        maxWavelengthIndex = np.where(abs(self.wavelengths-maxWavelength) ==
                               min(abs(self.wavelengths-maxWavelength)))[0]

        self.wavelengths = self.wavelengths[minWavelengthIndex:maxWavelengthIndex]
        self.spec = self.spec[minWavelengthIndex:maxWavelengthIndex]


        self.getSmoothedSpectrum(51)

    #LyAguess is in angstroms
    #check to see if the Lyman Alpha line is double peaked
    def isDoublePeaked(self, wavelengths, spec, LyAguess, noise):
        dspec = np.array([])
        dlambda = np.array([])

        #assume double peaked at first
        doublePeaked = True

        #set limits on the wavelength that we are looking at
        minWavelength = LyAguess-15
        maxWavelength = LyAguess+8
        minWavelengthIndex = np.where(abs(wavelengths-minWavelength) ==
                               min(abs(wavelengths-minWavelength)))[0][0]
        maxWavelengthIndex = np.where(abs(wavelengths-maxWavelength) ==
                               min(abs(wavelengths-maxWavelength)))[0][0]

        #we want to pick out just the lyman alpha line
        peakWavelengths = wavelengths[minWavelengthIndex:maxWavelengthIndex]
        peakSpec = spec[minWavelengthIndex:maxWavelengthIndex]
        maxNorm = max(peakSpec)
        peakSpec = peakSpec/maxNorm

        #take the derivative of the lyman alpha line profile
        for ii in range(len(peakWavelengths)-1):
            dspec = np.append(dspec, peakSpec[ii+1]-peakSpec[ii])
            dlambda = np.append(dlambda, (peakWavelengths[ii+1]+peakWavelengths[ii])/2.)
        
        #find which points have a negative slope
        negs=[]
        for jj in range(len(dspec)):
            if dspec[jj]<0:
                negs.append(jj)

        #find which points have a positive slope
        pos = [x+1 for x in list(set(range(len(spec)-1))-set(negs))]

        #points with a positive slope are shifted by one data point, 
        # then checked against negative slope points
        #common points will be the peaks
        peaks = list(set(pos) & set(negs))



        print("Npeaks: ", len(list(np.where(peakSpec[peaks]>0.8)[0])))
        print("Calculated peaks at: ", peakWavelengths[peaks])
        
        #first check for a large LyA peak that has a trough in it
        #this requires that for a spectrum normalized to the highest peak, 
        # the other part of the double peak(probably the blue peak) must
        # be at least 0.8.
        LyApeaks = np.array([])
        print("Number of peaks in the LyA is aaaa:{} higher: {}".format(len(peaks),
               len(np.where( peakSpec[peaks]>0.8 )[0])))
        print("List of peaks aaaa:{}".format(peakSpec[peaks]))
        print("List of high peaks aaaa:{}".format(np.where(peakSpec[peaks]>0.8)))
        if len(list(  np.where(  peakSpec[peaks]>0.8  )[0]  )) >= 2:
            for peak in list(np.where(peakSpec[peaks]>0.8)[0]):
                peak = peaks[peak]
                LyApeaks= np.append( LyApeaks,(np.where(wavelengths ==
                                       peakWavelengths[peak])[0][0]))
            #get an array of the indices of a sorted list
            indexSort = np.array(LyApeaks).argsort()
            #find the two highest peaks
            maxs = [indexSort[-2], indexSort[-1]]
            #find the index of the trough between the two max peaks
            trough = (self.spec[min(LyApeaks[maxs]):max(LyApeaks[maxs])].argmin()+
                      min(LyApeaks[maxs]))
            #check if the trough is significantly below the peaks
            #the trough must be less than "throughThresh" times the smaller of the peaks
            troughThresh = 0.8
            if self.spec[int(trough)] > troughThresh*min(self.spec[LyApeaks[maxs].astype(int)]):
                doublePeaked = False


            if doublePeaked:
                return trough 
        #only choose lines that are above the normalized noise level
        #this requires that if there is a smaller secondary peak, it must
        # be above the noise level(currently set at 2sigma), and the 
        # trough must be below that noise level
        LyApeaks=np.array([])
        nNoise = 3
        print("Checking for peaks above a threshold of {}".format(2*noise/maxNorm))
        if len(list(  np.where(  peakSpec[peaks]>  nNoise*noise/maxNorm  )[0]  )) >= 2:
            for peak in list(np.where(peakSpec[peaks]>nNoise*noise/maxNorm)[0]):
                peak = peaks[peak]
                print("Trying ", peakWavelengths[peak])
                print("Found peaks at: ", wavelengths[np.where(wavelengths ==
                                           peakWavelengths[peak])[0][0]])
                LyApeaks = np.append(LyApeaks, (np.where(wavelengths == 
                                           peakWavelengths[peak])[0][0]))
            #get an array of the indices of a sorted list
            indexSort = np.array(LyApeaks).argsort()
            #find the two highest peaks
            maxs = [indexSort[-2], indexSort[-1]]
            print(LyApeaks)
            print(min(LyApeaks[maxs]),max(LyApeaks[maxs]))
            #find the index of the trough between the two max peaks
            trough = (self.spec[min(LyApeaks[maxs]):max(LyApeaks[maxs])].argmin()+
                      min(LyApeaks[maxs]))

            #check if the trough is significantly below the peaks
            #the trough must be less than "throughThresh" which is the upper limit
            # for finding the secondary peaks 
            #some of these are divided by maxNorm because the noise level is already 
            # normalized, so I have to normalize the other comparisons
            troughThreshMin =2*noise/maxNorm 
            troughThreshMax = 0.5*min(  self.spec[LyApeaks[maxs].astype(int)]  )/maxNorm
            print("trough: {}, must be less than: {} and {}".format(self.spec[trough]/maxNorm,
                     troughThreshMax, troughThreshMin))
            if self.spec[trough]/maxNorm > troughThreshMin and \
                 self.spec[trough]/maxNorm > troughThreshMax:
                doublePeaked = False


            if doublePeaked:
                return trough
        return -2.0
    #get a smoothed spectrum
    def getSmoothedSpectrum(self, smoothFac):
        interp = interpolate.splrep(smoothSpec.smooth(self.wavelengths, smoothFac),
                                     smoothSpec.smooth(self.spec, smoothFac))
        lambdaNew = np.linspace(min(self.wavelengths), max(self.wavelengths), len(self.spec))
        self.smoothSpec = interpolate.splev(lambdaNew, interp, der=0)


    #remove a section of a spectrum
    #useful for atmospheric lines
    def removeSpec(self, center):
        print("Removing atmosphere line")
        #find index of dichroic boundary
        dichroIndex = np.where(abs(self.wavelengths-center)==
                              min(abs(self.wavelengths-center)))[0]
        #remove this boundary from the spectrum
        maskSize = 30
        for index in np.arange(dichroIndex-maskSize,dichroIndex+maskSize, 1):
            self.spec[index] = (self.spec[dichroIndex]+self.spec[dichroIndex+2*maskSize])/2.


    #find the rms error as a function of the spectrum
    #will be used to find where the emission lines are located in the spectra
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
        for ii in range(Nslice):
            rms = np.std(self.spec[ii*sliceSize:(ii+1)*sliceSize] -
                         self.smoothSpec[ii*sliceSize:(ii+1)*sliceSize])
            rmsarr = np.append(rmsarr, rms)
            self.peaks = np.append(self.peaks, (ii*sliceSize)+
                            np.where(self.spec[ii*sliceSize:(ii+1)*sliceSize]- 
                            self.smoothSpec[ii*sliceSize:(ii+1)*sliceSize] > 
                            rmsThresh*rms)[0])

            self.threshSpec = np.append(self.threshSpec, 
                        self.smoothSpec[ii*sliceSize:(ii+1)*sliceSize+1]+rmsThresh*rmsarr[ii])

        self.threshSpec = self.threshSpec[0:-1]

        self.emPeaks = (self.emRedshift+1)*LyA
        self.absPeaks = (self.absRedshift+1)*np.array([1303.3,1334.5])
        print("Found peaks at ", self.peaks)


    #if an emission line is large enough, it will contain more than 
    # one datapoint above the threshold for detecting lines
    #these lines will be counted more than once, so we remove all but the central data point.
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
            #if this point is in the same spectral line as the last point
            elif ( (self.wavelengths[self.peaks[jj]]-lastpeak) < 10):
                lastpeak = self.wavelengths[self.peaks[jj]]
                nInPeak += 1
            #if this point is in a new spectral line
            elif ( (self.wavelengths[self.peaks[jj]]-lastpeak) >= 10):
                lastpeak = self.wavelengths[self.peaks[jj]]
                #add the previous point to a list of line centers
                peakIndices = np.append(peakIndices, int((self.peaks[jj-1]+
                              self.peaks[jj-1-nInPeak])/2))
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


    #this will fit an emission line
    #returns the best fit center of the line
    def fitLine(self, peak, plotBool, subSpecN):
        #how many wavelength points away from center we want to consider
        #create subarrays from center +/- subSpecN points
        peak = peak[0]
        subSpectra = self.spec[peak-subSpecN:peak+subSpecN]
        subWavelengths = self.wavelengths[peak-subSpecN:peak+subSpecN]
        plotSpectra = self.spec[peak-3*subSpecN:peak+3*subSpecN]
        plotWavelengths = self.wavelengths[peak-3*subSpecN:peak+3*subSpecN]

        #calculate the noise of the spectrum
        noise = 1e-30*np.std(1e30*np.abs(self.spec[peak-2*subSpecN:peak-subSpecN]))
        #with subtracting out the continuum
        self.noiseSubtract =np.std(np.abs(self.spec[peak-2*subSpecN:peak-subSpecN]-
                                    self.smoothSpec[peak-2*subSpecN:peak-subSpecN])) 

        #check if the line is double peaked, then it will have different 
        # behavior when calculating systematic redshift
        trough = self.isDoublePeaked(self.wavelengths, self.spec, self.wavelengths[peak], noise)
        if trough > 0:
            doublePeaked = True
            self.peakType = 'd'    
        else:
            doublePeaked = False

#        if len(doublePeaks) >= 2:
#            doublePeaked = True
#        else:
#            doublePeaked = False
#
#        #the redshift will then be fit to the trough in the center
#        if doublePeaked:
#            #get an array of the indices of a sorted list
#            indexSort = doublePeaks.argsort()
#            #find the two highest peaks
#            maxs = [indexSort[-2], indexSort[-1]]
#            #find the index of the trough between the two max peaks
#            trough = self.spec[min(doublePeaks[maxs]):max(doublePeaks[maxs])].argmin()+min(doublePeaks[maxs])
#
#            #check if the trough is significantly below the peaks
#            #the trough must be less than "throughThresh" times the smaller of the peaks
#            troughThresh = 0.8
#            if self.spec[trough] > troughThresh*min(self.spec[doublePeaks[maxs]]):
#                doublePeaked = False

        #number of elements in the spectrum array
        N = np.shape(self.spec)[0]
        #try to fit to gaussian
        #first, check to see if the line is close to the edge
        if (peak > subSpecN) and (peak < N-subSpecN):
            maxSpec = max(subSpectra)

            #perform fit of gaussian
            fitParam, fitCov = curve_fit(gauss, (subWavelengths-
                                subWavelengths[subSpecN-1]), subSpectra/maxSpec)
            #print out fit data
            print("A = ", fitParam[0])
            print("Center = ", fitParam[1]+subWavelengths[subSpecN-1], "Angstroms")
            print("Sigma = ", fitParam[2])

            #set all of the line fit parameters to member variables
            self.lineAmplitude = fitParam[0]
            self.lineCenter = fitParam[1]+subWavelengths[subSpecN-1] #in angstroms
            self.lineSigma = fitParam[2]
            self.lineConst = fitParam[3]
            self.lineLinear = fitParam[4]

            #calculate redshift
            z = (fitParam[1]+subWavelengths[subSpecN-1])/LyA - 1

            #create a plot of each spectrum and the fitted line.
            if (plotBool):
                #putting this set of images in a folder for each object
                print("Checking to see if directory exists")
                if not(os.path.isdir("./images/"+self.objName)):
                    os.makedirs("./images/"+self.objName)
                #check if the feature is emission
                if fitParam[0]>0:
                    #comment out to save figures sorted by object
                    plt.cla()
                    #plot the spectrum
                    plt.plot(self.wavelengths, self.spec)
                    #save a figure of the entire spectrum
                    plt.savefig("./images/"+self.objName+"/"+self.filename[0:-5]+"fullspec.png")
                    

                    plt.cla()
                    #plot the data near the LyA line
                    plt.plot(subWavelengths, subSpectra/maxSpec, linewidth=2)
                    #plot the best fit
                    plt.plot(subWavelengths, (subWavelengths-subWavelengths[subSpecN-1])*fitParam[4]+fitParam[3]+fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[subSpecN-1]))**2/(2*fitParam[2]**2)), linewidth=2, label="z=%.3f" % z+"d"*doublePeaked)
                    plt.legend()
                    #save the figure of just the LyA line
                    plt.savefig("./images/"+self.objName+"/"+self.filename[0:-5]+"linespec.png")

                    #now save a copy of the images all together in one directory, prepended by the redshift
                    if not(os.path.isdir("./images/all/")):
                        os.makedirs("./images/all/")
                    plt.cla()
                    #plot the spectrum
                    plt.plot(self.wavelengths, self.spec)
                    #save the figure of the whole spectrum
                    plt.savefig("./images/all/{:6.4f}_".format(z)+self.objName+"_"+self.filename[0:-5]+"_fullspec.png")
                    
                    plt.cla()
                    #plot the data around LyA
                    plt.plot(subWavelengths, subSpectra/maxSpec, linewidth=2)
                    #plot the fit of LyA
                    plt.plot(subWavelengths, (subWavelengths-subWavelengths[subSpecN-1])*fitParam[4]+fitParam[3]+fitParam[0]*np.exp(-(subWavelengths-(fitParam[1]+subWavelengths[subSpecN-1]))**2/(2*fitParam[2]**2)), linewidth=2, label="z=%.3f" % z+"d"*doublePeaked)


                    plt.legend()
                    plt.xlim([subWavelengths[0], subWavelengths[-1]])



                    #if it is a double peak, we will plot additional things
                    if doublePeaked:
                        #plot the location of the trough
                        plt.plot([self.wavelengths[trough], self.wavelengths[trough]], [0,1], 'k:', linewidth=2)
                    plt.savefig("./images/all/{:6.4f}".format(z)+
                            "d"*doublePeaked+"_"+self.objName+"_"+self.filename[0:-5]+"_linespec.png")


                        #plot the position of the peaks in the LyA line
#                        plt.plot(self.wavelengths[doublePeaks], self.spec[doublePeaks]/max(subSpectra), 'ko')


                    #save the figure with just the LyA line
#                    plt.savefig("./images/all/{:6.4f}".format(z)+
#                        "d"*doublePeaked+"_"+self.filename[0:-5]+"_linespec.png")

            #if it is not double peaked, we will find the equivalent width
            if not doublePeaked:
                #This sets the width of the line to 3 sigma
                width = 3*fitParam[2]
                #this finds the center of the line
                center = fitParam[1]+subWavelengths[subSpecN-1]
                #calculate an array of wavelengths at the center +/- width
                continuumWavelengths = subWavelengths.clip((center-width), (center+width))
                #get the data in the same range
                continuumSpec = fitParam[3]+fitParam[4]*(continuumWavelengths-
                                    subWavelengths[subSpecN-1])
                #this is the fit to the continuum on the line
                linSpec = self.spec[np.where(subWavelengths.clip(center-width, center+width)==
                            continuumWavelengths)]

                #find the equivalent width of LyA if not double peaked
                self.EW = np.trapz((continuumSpec-linSpec)/continuumSpec, x=continuumWavelengths)
                print("Calculated the equivalent width to be: {:5.3f} Angstroms".format(self.EW))
                # if fitParam[0]<0:
                #     return -1*(fitParam[1]+subWavelengths[subSpecN-1])/LyA+1

                #make sure this is not fitting any absorption noise features
                if fitParam[0]>0:
                    return (fitParam[1]+subWavelengths[subSpecN-1])/LyA-1
                else:
                    return -3
            else:
                #if the line id double peaked
                return self.wavelengths[trough]/LyA - 1


        else:
            return -2


    #this will fit a series of absorption lines, requires a guess of the redshift
    #uses cross correlation to find the most likely absorption redshift
    def fitAbsLine(self, zguess, plotBool):
        print("Calculating absorption redshift")
        #zguess = self.wavelengths[peak[0]]/LyA-1
        print("Zguess is ", zguess)

        #this is how far away from a guess we want to check
        absZthresh = 0.005

        #create a list of redshifts to iterate over
        zlist = np.linspace(zguess*(1-absZthresh),zguess*(1+absZthresh), 100)
        integrals = np.array([])
        zstep = 0

        #for each iterative redshift
        for z in zlist:
            lineSpec = np.zeros(np.shape(self.fullspec)[0])

            #for each metal line in our list, create a line and add it to our comparison spectrum
            for kk in range(lines.size):
                line = ((z+1)*lines[kk])
                gaussian = f[kk]*np.exp(-((self.fullwavelengths - line)**2)/(2*W[kk]**2))
                lineSpec += gaussian

            #multiply our comparison spectrum by the data, then integrate
            multSpec = (self.fullspec)*lineSpec
            integrals = np.append(integrals, np.trapz(self.fullwavelengths, multSpec))

        #pick out the peak of the cross correlation
        absZ = zlist[integrals.argmax()]

        #if wanted, create a plot of the spectrum with absorption lines noted.
        if (plotBool):
            print("Checking to see if directory exists")
            if not(os.path.isdir("./images/"+self.objName)):
                os.makedirs("./images/"+self.objName)
                print("Created directory")

            plt.cla()
            plt.plot(self.fullwavelengths, self.fullspec/max(self.spec))
            #plot each of the absorption lines
            for ii in range(len(lines*(absZ+1))):
                peak = (lines*(absZ+1))[ii]
                plt.plot([peak, peak], [-1, 1], 'r--', linewidth = 2)


            plt.xlim([min(lines*(absZ+1)), max(lines*(absZ+1))])
            plt.ylim([-1, 1])
            plt.savefig("./images/"+self.objName+"/"+self.filename[0:-5]+"absorbspec.png")

            if not(os.path.isdir("./images/all/absorption/")):
                os.makedirs("./images/all/absorption/")
            plt.cla()
            plt.plot(self.fullwavelengths, self.fullspec/max(self.spec))
            #plot each of the absorptoin lines
            for ii in range(len(lines*(absZ+1))):
                peak = (lines*(absZ+1))[ii]
                plt.plot([peak, peak], [-1, 1], 'r--', linewidth = 2)


            plt.xlim([min(lines*(absZ+1)), max(lines*(absZ+1))])
            plt.ylim([-1, 1])
            plt.savefig("./images/all/absorption/"+str(absZ)+"_"+self.objName+"_"+self.filename[0:-5]+"_absorbspec.png")
        if (absZ < 3.35) and (absZ > 2.7):
            return absZ
        else:
            return -2


    #find the equivalent width of the Lyman alpha line
    def calculateEQ(self):
        return
