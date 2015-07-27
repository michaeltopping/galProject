#Will contain cosmology functions such as distance conversions.
#Distances based on "Distance measures in cosmology" David Hogg 2000
import numpy as np
import scipy

#Angular size distance
def DA(z):
	return DM(z)/(1+z)

#Comoving distance (transverse)
def DM(z, OmegaM=0.308, OmegaL=0.692):
	OmegaK = 1-OmegaM-OmegaL
	if (OmegaK>0):
		return DH(z)/np.sqrt(Omegak(z))*np.sinh(np.sqrt(OmegaK(z))*DC(z)/DH(z)) 
	elif (OmegaK == 0):
		return DC(z)
	elif (OmegaK < 0):
		return DH(z)/np.sqrt(-Omegak(z))*np.sin(np.sqrt(-OmegaK(z))*DC(z)/DH(z)) 
	else:
		raise ValueError('Invalid value for OmegaK')


#Hubble Distance
def DH(H0=67.8):
	c = 3e10
	#unit conversion of hubble constant
	H0 = H0*3.214e-20
	return c/H0


#Coming distance (line-of-sight)
def DC(z, OmegaM=0.308, OmegaL = 0.692):
	OmegaK = 1-OmegaM-OmegaL
	Ez = lambda z: 1/np.sqrt(OmegaM*(1+z)**3+OmegaK*(1+z)**2+OmegaL)
	return DH()*scipy.integrate.quad(Ez, 0, z)[0]
