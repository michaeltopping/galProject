import numpy as np

filename = "../spec/notes/all.cat"

data = np.loadtxt(filename, dtype="string")

print np.shape(data)