import numpy as np 
import os

objects = np.array([])

#find all filenames in the directory 
for filename in os.listdir("../spec/mask_design/"):
	#we only want the files that are used to generate slitmasks
	if filename.endswith(".autoslitin"):
		#read in the data
		data = np.genfromtxt("../spec/mask_design/"+filename, skip_header=1, dtype="str")
		#for each row in the data
		for ii in range(np.shape(data)[0]):
			#sed the object name
			obj = data[ii][0]
			#set the RA and dec in degrees
			RA = float(data[ii][3])*15+float(data[ii][4])*15./60. + float(data[ii][5])*15/60./60.
			dec = float(data[ii][6])+float(data[ii][7])/60.+float(data[ii][8])/60./60.
			if objects.size == 0:
				newline = [obj, RA, dec]
				objects = np.reshape(np.append(objects, newline), (-1,3))
			elif not (obj in objects[:,0]):
				newline = [obj, RA, dec]
				objects = np.reshape(np.append(objects, newline), (-1,3))
			else:
				print("Found Duplicate: {}".format(obj))



#now add in objects from other sources
data = np.genfromtxt("./positions_H.dat", dtype="str")
for ii in range(np.shape(data)[0]):
	obj = data[ii][0]
	RA_H = data[ii][1]
	dec_h = data[ii][2]
	RA = (int(RA_H[0:2])+int(RA_H[3:5])/60.+float(RA_H[6:13])/3600.)*15
	dec = (int(dec_h[0:1])+int(dec_h[2:4])/60.+float(dec_h[5:13])/3600.)

	if not (obj in objects[:,0]):
		newline = [obj, RA, dec]
		objects = np.reshape(np.append(objects, newline), (-1,3))
	else:
		print("Found Duplicate: {}".format(obj))

if __name__ == "__main__":
	#print(objects)
	masterFile = open("masterPositions.dat", "w")
	for row in objects[:]:
		print row
		masterFile.write("{:11s}    {:9.5f}    {:9.5f}\n".format(row[0], float(row[1]), float(row[2])))
	masterFile.close()
