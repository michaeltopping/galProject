import numpy as np
from math import *

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
		self.RA = 0
		self.dec = 0


	#add a redshift measurement
	#type deciedes whether it is from emission lines, or absorption lines
	def addRedshift(self, type, z):
		if type=="em":
			self.emRedshifts = np.append(self.emRedshifts, z)
		elif type=="abs":
			self.absRedshifts = np.append(self.absRedshifts, z)
		else:
			print("Incorrect type, no redshifts was added.")

	#add a type to the galaxy
	#these types include lae, lbg, double peaked...
	def addType(self, type):
		self.type += type

	#change the redshift based on the type of spectral features it contains.
	# shift data from Adelberger et al. 2003
	def systematicShift(self):

		#if the lyman alpha line is double peaked
		if ("d" in self.type):
			z = np.average(self.emRedshifts)
			self.sysRedshift = z

		#if the galaxy has only LyA emission, and no absorption features.
		#we want to shift it by 310km/s on average
		elif (self.emRedshifts.size > 0 and self.absRedshifts.size == 0):
			z = np.average(self.emRedshifts)
			#shift redshift into velocity space, then add the 310km/s then shift back to redshift space
			vel = ((z+1)**2-1)/((z+1)**2+1)
			vel = (vel*3e5-310.)/3e5
			z = sqrt((1+vel)/(1-vel))-1
			self.sysRedshift = z


		#if the galaxy shows LyA emission features as well as interstellar absorption features.
		elif (self.emRedshifts.size > 0 and self.absRedshifts.size > 0):
			zEm = np.average(self.emRedshifts)
			zAbs = np.average(self.absRedshifts)
			velEm = ((zEm+1)**2-1)/((zEm+1)**2+1)
			velAbs = ((zAbs+1)**2-1)/((zAbs+1)**2+1)
			vel = (velEm*3e5-230.+0.114*abs(velEm-velAbs))/3e5
			z = sqrt((1+vel)/(1-vel))-1
			self.sysRedshift = z

		#if the galaxy only shows interstellar absorption features.	
		elif (self.absRedshifts.size > 0 and self.emRedshifts.size == 0):
			z = np.average(self.absRedshifts)
			vel = ((z+1)**2-1)/((z+1)**2+1)
			vel = (vel*3e5+150.)/3e5
			z = sqrt((1+vel)/(1-vel))-1
			self.sysRedshift = z

		#if for some reason things don't work, this will happend
		else:
			print("Error calculating systematic shift, no lines detected.")
			self.sysRedshift = np.average(self.emRedshifts)



	def setCoords(self, RA, dec):
		self.RA = RA
		self.dec = dec
