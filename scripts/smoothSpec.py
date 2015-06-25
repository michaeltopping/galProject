#Smooth out an array
import numpy as np

#input:
# array: array that you want smoothed
# factor: size of bins that you want for the smoothing.  Must be an odd number
def smooth(array, factor):
	if (factor%2) == 0:
		raise(ValueError('Must smooth by an odd number'))

	newArr = np.array([])
	for ii in range(int(len(array)/factor)):
		newArr = np.append(newArr, np.average(array[ii*factor:ii*factor+factor]))
	return newArr


#for testing
if __name__=="__main__":
	print(smooth([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 3))

	
