import numpy as np
import random
from math import *


#constants
c = 3e10

class Photon(object):
    """docstring for the Photon class"""
    def __init__(self, nu):
        super(, self).__init__()
        self.nu = nu
        self.theta = 2*pi*random.random()
        self.phi = 2*pi*random.ramdom()

    def dopplerShift(self, vel):
        self.nu *= vel/c  #fix this
