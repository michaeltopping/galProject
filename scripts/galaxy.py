import numpy as np

#constants
LyA = 1215.67

#galaxy class
#will contain the redshift data for multiple observations
class galaxy():
	def __init__(self):
		self.emRedshifts = np.array([])
		self.absRedshifts = np.array([])
		self.type = ""
		self.sysRedshift = 0


	def addRedshift(self, type, z):
		if type=="em":
			self.emRedshifts = np.append(self.emRedshifts, z)
		elif type=="abs":
			self.absRedshifts = np.append(self.absRedshifts, z)
		else:
			print("Incorrect type, no redshifts was added.")

	def addType(self, type):
		self.type += type

	def systematicShift(self):
		#change the redshift based on the type of spectral features it contains.
		# shift data from Adelberger et al. 2003
		if (self.emRedshifts.size > 0 and self.absRedshifts.size == 0):
			self.sysRedshift = np.average(self.emRedshifts)
		elif (self.emRedshifts.size > 0 and self.absRedshifts.size > 0):
			self.sysRedshift = (np.average(self.emRedshifts) + self.absRedshifts)/2.
		elif (self.absRedshifts.size > 0 and self.emRedshifts.size == 0):
			self.sysRedshift = self.absRedshifts
