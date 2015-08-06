import numpy as np
import matplotlib.pyplot as plt

# fileout = open("nbpositions.dat", 'w')
# filenames = ['ssa22a_Lnb1.autoslitin','ssa22a_Lnb2.autoslitin']
# for filename in filenames:
# 	data = np.genfromtxt("../spec/mask_design/"+filename, skip_header=1, dtype="str")
# 	#for each row in the data
# 	for ii in range(np.shape(data)[0]):
# 		#sed the object name
# 		obj = data[ii][0]
# 		#set the RA and dec in degrees
# 		fileout.write("{}    {}:{}:{}    {}:{}:{}\n".format(obj,   data[ii][3],data[ii][4],data[ii][5],data[ii][6],data[ii][7],data[ii][8],))

# fileout.close()


#read in the data from the table of objects
data = np.genfromtxt("./nbpositions_pixel.dat", dtype="|S9")

#set the columns to appropriate arrays
names = data[:,0].astype(str)
RA = data[:,1].astype(float)
dec = data[:,2].astype(float)
for ii in range(len(names)):
	print("{} {} {}".format(names[ii], RA[ii], dec[ii]))
centroid = True


if centroid == True:
	#read in imcentroid output
	data = np.genfromtxt("./nbpositions_centroid.dat")
	iis = data[:,5]
	xs = data[:,1]
	ys = data[:,3]
	regionFile = open("nbpositions_adjusted.dat", "w")
	difx = []
	dify = []

#find the average shift
# The following is an if statement
if centroid == True:
	# The following is a for loop
	for ii in range(len(names)):
		# The following is an if statement
		if ii in iis:
			index = np.where(iis == ii)

			difx.append(float(xs[index])-RA[ii-1])
			dify.append(float(ys[index])-dec[ii-1])
	dx = np.average(difx)
	dy = np.average(dify)
# Now we gonna print some shit
print("Average shift: dx={}  dy={}".format(dx, dy))
# Woohoo another for loop
for ii in range(len(names)+1):

	if centroid == True:
		if ii in iis:
			index = np.where(iis == ii)
			print("{} {}".format(xs[index], RA[ii-1]))

			regionFile.write("{} {} {} \n".format(names[index][0], float(xs[index]), float(ys[index])))
		else:
			regionFile.write("{} {} {} \n".format(names[ii-1], RA[ii-1]+dx, dec[ii-1]+dy))
	elif centroid == False:
		regionFile.write("{} {} {} \n".format(RA[ii], dec[ii]))

print(difx)

