#Testing detecting spectral lines using a neural network
#Test spectra objects:
#  ssaly_3.b.NB407.msdfc_v
#  ssaly_3.b.NB2398.msdfc_v
#
#
#




import numpy as np 
import pyfits
import sys
import time

#definition of sigmoid function, optional return of derivative
def sigmoid(x, deriv=False):
	#this will return the derivative at x
	if deriv:
		return x*(1-x)
	#otherwise, return the value of the function at x
	return 1/(1+np.exp(-x))



spec, header = pyfits.getdata("../spec/ssaly_3.b.NB407.msdfc_v.fits", 0, header=True)

#creating wavelength solution
lambdaMin = float(header['CRVAL1'])
dlambda = float(header['CD1_1'])
Nlambda = float(header['NAXIS1'])

wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)



#set limits on the wavelength that we are looking at
minWavelength = 4500
maxWavelength = 5300
minWavelengthIndex = np.where(abs(wavelengths-minWavelength) == min(abs(wavelengths-minWavelength)))[0]
maxWavelengthIndex = np.where(abs(wavelengths-maxWavelength) == min(abs(wavelengths-maxWavelength)))[0]

wavelengths = wavelengths[minWavelengthIndex:maxWavelengthIndex]
spec = spec[minWavelengthIndex:maxWavelengthIndex]



nTrain = 1
specSize = spec.size
#define the matrix of training data
X = np.array([spec])

#define the vector of solutions
Y = np.array([1])

#seed the random number generators
np.random.seed(6563)

#randomly initialize the synapse weight matrix
#were doing a 3-layer neural network so that we can 
#detect correlations between different data points
synapse0 = 2*np.random.random((specSize,nTrain)) - 1
synapse1 = 2*np.random.random((nTrain,1)) - 1

#begin training the weight matrices
for ii in xrange(60000):

	#progress bar
	sys.stdout.write("\rTraining: {}%".format((ii+1)/600))
	sys.stdout.flush()

	#put the data through the layers of the neural network.
	layer0 = X
	layer1 = sigmoid(np.dot(layer0,synapse0))
	layer2 = sigmoid(np.dot(layer1,synapse1))

	#calculate the error from the weights
	layer2Error = Y - layer2

	#find in what direction our guess was off
	layer2Delta = layer2Error*sigmoid(layer2,deriv=True)

	#calculate how much of this error is attributed to layer1 on layer2
	layer1Error = layer2Delta.dot(synapse1.T)

	#now find which direction the layer1 errors are off from
	layer1Delta = layer1Error * sigmoid(layer1, deriv=True)

	#update the weights
	synapse1 += layer1.T.dot(layer2Delta)
	synapse0 += layer0.T.dot(layer1Delta)


print("--------Test Cases--------")
spec, header = pyfits.getdata("../spec/ssaly_3.b.NB2398.msdfc_v.fits", 0, header=True)

#creating wavelength solution
lambdaMin = float(header['CRVAL1'])
dlambda = float(header['CD1_1'])
Nlambda = float(header['NAXIS1'])

wavelengths = np.arange(lambdaMin, lambdaMin + int(dlambda*(Nlambda)), dlambda)



#set limits on the wavelength that we are looking at
minWavelength = 4500
maxWavelength = 5300
minWavelengthIndex = np.where(abs(wavelengths-minWavelength) == min(abs(wavelengths-minWavelength)))[0]
maxWavelengthIndex = np.where(abs(wavelengths-maxWavelength) == min(abs(wavelengths-maxWavelength)))[0]

wavelengths = wavelengths[minWavelengthIndex:maxWavelengthIndex]
spec = spec[minWavelengthIndex:maxWavelengthIndex]
print("Likelyhood of line existence:")
x1 = np.array(spec)
print("Spectra: {}".format(sigmoid(np.dot(np.dot(x1.T,synapse0),synapse1))))
x1 = np.random.rand(spec.size)
print("Spectra: {}".format(sigmoid(np.dot(np.dot(x1.T,synapse0),synapse1))))

