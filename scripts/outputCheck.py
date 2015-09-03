import numpy as np




#define the data that is output from analysis
analysisOutput = np.genfromtxt("./table.dat", dtype=str, skip_header=1)
#define dictionary of objects output from analysis
analysisDict = {}
analysisDictAbs = {}

#create the data from the guesses
inputOldmask = np.genfromtxt("./ssa22a.oldmaskkey.txt", dtype=str)
inputCurrent = np.genfromtxt("../spec/mask_design/z.ssa22a.current", dtype=str, 
                                skip_header=1)
inputOrig = np.genfromtxt("../spec/all.cat", dtype=str)
inputDict = {}
inputDictAbs = {}


#loop through each object in the analysis output
for objName, zem, zabs in zip(analysisOutput[:,0], analysisOutput[:,2], analysisOutput[:,3]):
    analysisDict[objName] = float(zem)
    analysisDictAbs[objName] = float(zabs)


#loop through each object in the input files
for objName, zem, zabs in zip(inputOldmask[:,0], inputOldmask[:,5], inputOldmask[:,6]):
    inputDict[objName[7:]] = float(zem)
    inputDictAbs[objName[7:]] = float(zabs)

for objName, zem, zabs in zip(inputCurrent[:,2], inputCurrent[:,8], inputCurrent[:,9]):
    inputDict[objName] = float(zem)
    inputDictAbs[objName] = float(zabs)

for objName, zem, zabs in zip(inputOrig[:,0], inputOrig[:,2], inputOrig[:,3]):
    inputDict[objName] = float(zem)
    inputDictAbs[objName] = float(zabs)

print("Begin checking {} input objects against {} measured objects.".format(len(inputDict), len(analysisDict)))
print("----Emission Redshifts----")

#loop through each object in the input and print out the measured value
print("Objects that do not have measured redshifts")
print("object    guess    measured")
for obj in sorted(inputDict):
    if inputDict[obj] > 0:
        if not obj in analysisDict:
            continue
        if analysisDict[obj] < 0:
            print("{}    {}    {}".format(obj, inputDict[obj], analysisDict[obj]))

print()
print("Objects that have differing redshifts")
threshold = 0.003
print("object    guess    measured")
for obj in sorted(inputDict):
    if inputDict[obj] > 0:
        if not obj in analysisDict:
            analysisDict[obj] = -7.0
        if analysisDict[obj] > 0:
            if abs(analysisDict[obj] - inputDict[obj])>threshold:
                print("{}    {}    {}".format(obj, inputDict[obj], analysisDict[obj]))



print("----Absorption Redshifts----")
#loop through each object in the input and print out the measured value
print("Objects that do not have measured redshifts")
print("object    guess    measured")
for obj in sorted(inputDictAbs):
    if inputDictAbs[obj] > 0:
        if not obj in analysisDictAbs:
            continue
        if inputDictAbs[obj] > 2.3: 
            if analysisDictAbs[obj] < 0:
                print("{}    {}    {}".format(obj, inputDictAbs[obj], analysisDictAbs[obj]))

print()
print("Objects that have differing redshifts")
threshold = 0.003
print("object    guess    measured")
for obj in sorted(inputDictAbs):
    if inputDictAbs[obj] > 0:
        if not obj in analysisDictAbs:
            analysisDictAbs[obj] = -7.0
        if analysisDictAbs[obj] > 0:
            if inputDictAbs[obj] > 2.3:
                if abs(analysisDictAbs[obj] - inputDictAbs[obj])>threshold:
                    print("{}    {}    {}".format(obj, inputDictAbs[obj], analysisDictAbs[obj]))





#now print a list of all objects
#
#for obj in sorted(analysisDict):
#    if obj in inputDict:
#        if inputDict[obj] > 0:
#            print("{}  {}  {}   {}".format(obj, inputDict[obj], analysisDict[obj], analysisDict[obj]-inputDict[obj]))
#


