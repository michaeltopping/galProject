import numpy as np
from math import *

#constants
LyA = 1215.67

#object class
#used as a member of dataset for tracking which objects belong with which dataset
class Object():
	def __init__(self):
		self.emRedshifts = np.array([])
		self.absRedshifts = np.array([])
		self.type = ""
		self.sysRedshift = 0
		self.RA = 0
		self.dec = 0


class Galaxy():
	def __init__(self, name):
		self.name = name
		self.emRedshifts = []
		self.absRedshifts = []
		self.type = []
		self.sysRedshift = -9.9999


	#change the redshift based on the type of spectral features it contains.
	# shift data from Adelberger et al. 2003, Trainor et al. 2015
	def systemicShift(self):

		#if the lyman alpha line is double peaked
		if ("d" in self.type):
			z = np.average(self.emRedshifts)
			self.sysRedshift = z

		#if the galaxy has only LyA emission, and no absorption features.
		#we want to shift it by 200m/s on average
		elif (len(self.emRedshifts) > 0 and len(self.absRedshifts) == 0):
			if self.name.startswith("C") or self.name.startswith("D") or self.name.startswith("M") or self.name.startswith("MD"):
				z = np.average(self.emRedshifts)
				z -= (310./3e5)*(1+z)
			else:
				z = np.average(self.emRedshifts)
				z -= (200./3e5)*(1+z)
			self.sysRedshift = z


		#if the galaxy shows LyA emission features as well as interstellar absorption features.
		elif (len(self.emRedshifts) > 0 and len(self.absRedshifts) > 0):
			zEm = np.average(self.emRedshifts)
			zAbs = np.average(self.absRedshifts)
			dv = (zEm-zAbs)/zEm*3e5
			zEm -= (230.-0.114*dv)/3e5*(1+zEm)
			self.sysRedshift = zEm

		#if the galaxy only shows interstellar absorption features.	
		elif (len(self.absRedshifts) > 0 and len(self.emRedshifts) == 0):
			z = np.average(self.absRedshifts)
			z += (150./3e5)*(1+z)
			self.sysRedshift = z

		#if for some reason things don't work, this will happend
		else:
			self.sysRedshift = -9.9999

		self.emRedshift = np.average(self.emRedshifts)
		if np.isnan(self.emRedshift):
			self.emRedshift = -9.9999
		self.absRedshift = np.average(self.absRedshifts)
		if np.isnan(self.absRedshift):
			self.absRedshift = -9.9999
