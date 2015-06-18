import numpy as np


#galaxy class
#will contain the redshift data for multiple observations
class galaxy():
	def __init__(self):
		self.emRedshifts = np.array([])
		self.absRedshifts = np.array([])
		
		
	def addRedshift(self, type, z):
		if type=="em":
			self.emRedshifts = np.append(self.emRedshifts, z)
		elif type=="abs":
			self.absRedshifts = np.append(self.absRedshifts, z)
		else:
			print("Incorrect type, no redshifts was added.")
		