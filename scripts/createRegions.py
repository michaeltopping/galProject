import numpy as np
import matplotlib.pyplot as plt


#read in the data from the table of objects
data = np.genfromtxt("../spec/mask_design/z.ssa22a.current_pixel", dtype="|S9")

#set the columns to appropriate arrays
names = data[:,2].astype(str)
RA = data[:,11].astype(float)
dec = data[:,12].astype(float)
print(names)
centroid = True


if centroid == True:
	#read in imcentroid output
	data = np.genfromtxt("./centroid.dat")
	iis = data[:,5]
	xs = data[:,1]
	ys = data[:,3]
	regionFile = open("positions.dat", "w")
	difx = []
	dify = []

#find the average shift
if centroid == True:
	for ii in range(len(names)):
		if ii in iis:
			index = np.where(iis == ii)

			difx.append(float(xs[index])-RA[ii-1])
			dify.append(float(ys[index])-dec[ii-1])
	dx = np.average(difx)
	dy = np.average(dify)
print("Average shift: dx={}  dy={}".format(dx, dy))
for ii in range(len(names)+1):

	if centroid == True:
		if ii in iis:
			index = np.where(iis == ii)
			regionFile.write("{} {} {} \n".format(names[index][0], float(xs[index]), float(ys[index])))
		else:
			regionFile.write("{} {} {} \n".format(names[ii-1], RA[ii-1]+dx, dec[ii-1]+dy))
	elif centroid == False:
		regionFile.write("{} {} {} \n".format(RA[ii], dec[ii]))


plt.hist(difx, bins=30, histtype="stepfilled", color='blue', alpha=0.5, range=(-2, 7))
plt.hist(dify, bins=30, histtype="stepfilled", color='red', alpha=0.5, range=(-2, 7))
plt.show()

