import numpy as np
import operator
from math import *
import matplotlib.pyplot as plt
from cosmology import *



# create the list of Millennium simulation timesteps
redshiftList = open("./halos/redshift_text.txt", 'r')
redshifts = []
for r in redshiftList:
    redshifts.append(float(r))

# reverse the redshift list
redshifts = redshifts[::-1]

# physical constants
c = 3e10
sToYr = (1./3.15e7)



class halo():

    def __init__(self, row):
        # assign the variables
        self.haloid = row[2]
        self.z = float(row[14])
        self.desid = row[5]
        self.nextHaloId = row[9]
        self.firstProgId = row[10]
        self.nextProgId = row[11]
        self.m = 8.6e8*float(row[15])

    def set_des(self, halo):
        self.des = halo






class fof_group():
    
    def __init__(self, id):
        self.halos = []
        self.m = 0
        self.id = id


    def find_total_mass(self):
        for h in self.halos:
            self.m += h.m 


    # sort the halos by their mass
    def sort_halos(self):
       self.halos.sort(key=lambda x: x.m, reverse=True) 


    # fill the halos array with all halos in the fof group
#    def fill_halos(self, halos, firstProgId):
#        nextProgId = 0
#        if firstProgId in halos:
#            self.halos.append(halos[firstProgId])
#            
#            nextProgId = halos[firstProgId].nextProgId
#            del halos[firstProgId] 
#
#
#        if nextProgId in halos:
#            print("NextProgId is in halos")
#            lastProg = False
#        else:
#            print("The next progenitor id is not in halos")
#            lastProg = True
#
#        # fill the rest of the halos
#        while not lastProg:
#            self.halos.append(halos[nextProgId])
#
#            
#            # set the next progenitor id
#            oldnextProgId = nextProgId
#            nextProgId = halos[nextProgId].nextProgId
#            del halos[oldnextProgId]
#            if nextProgId in halos:
#                lastProg = False
#            else:
#                lastProg = True
#
#
#        return halos
#
#
#
#    def add_halo(self, halos, nextProgId):
#        
#        self.halos.append(halos[nextProgId])
#
#        nextProgId
        





class snapnum():
    def __init__(self, snapId):
        self.snapId = snapId
        self.fofGroups = {}



    def create_fof(self, haloList, firstFofId):
        # make sure the fof group is not repeated in the fofGroups dictionary
        #  the fof groups are labelled by their first group member ID
        if not firstFofId in self.fofGroups:

            self.fofGroups[firstFofId]=fof_group(firstFofId)

            # loop through the halolist for this redshift/snapnum
            #  and add to the group every halo that references the firstFOF
            #  as as a member of the FOF group.
            for haloKey in haloList:
                halo = haloList[haloKey]
                if halo.firstProgId == firstFofId:
                    # add it to the fof Groups list
                    self.fofGroups[firstFofId].halos.append(halo)

        self.fofGroups[firstFofId].find_total_mass()
        self.fofGroups[firstFofId].sort_halos()



    def get_nHalos(self):
        # get the number of halos in this snapshot
        nHalos = 0
        for key in self.fofGroups:
            nHalos += len(self.fofGroups[key].halos)

        # assign the value to a memeber variable
        self.nHalos = nHalos
        self.nGroups = len(self.fofGroups)



    def arrange_xs(self, maxgroups):
        ii = 0
        self.groupLocs = {}
        self.xs = np.arange((maxgroups-self.nGroups)/2., (maxgroups+self.nGroups)/2.)
        for fof in (sorted(self.fofGroups.values(), key=operator.attrgetter('m'), reverse=True)):
            self.groupLocs[fof.halos[0].haloid] = (self.xs[ii], self.snapId)
            ii += 1
#


    def plot_snapshot(self):
        ii = 0
        for fof in (sorted(self.fofGroups.values(), key=operator.attrgetter('m'), reverse=True)):
#            fof = self.fofGroups[key]
            plt.plot(self.xs[ii], DC(redshifts[self.snapId])*sToYr/c, 'ko', markersize = (max(sqrt(fof.m), 5e6))*5e-7)
            ii += 1



    def plot_prog(self, lastsnap):

        ii = 0
        # loop through each of the fof groups in this snapshot
        for fof in (sorted(self.fofGroups.values(), key=operator.attrgetter('m'), reverse=True)):
            nextProg = fof.halos[0].desid

            print("Looking for id: {} in snapshot: {}, from halo: {}".format(nextProg, self.snapId, fof.halos[0].haloid))
            # now loop through each of the halos in the next snapshot to find where the halo goes into
            
            for nextfof in (sorted(lastsnap.fofGroups.values(), key=operator.attrgetter('m'), reverse=True)):
                firstHaloId = nextfof.halos[0].haloid
                # loop through each of the halos in the fof group
                #  look for any halos in the group that will be desendants of previous halos
                for subhalo in nextfof.halos:
                    if int(subhalo.haloid) == int(nextProg):
                        nextProg = float("{:.1f}".format(firstHaloId))
                        print("Found match in snap: {} for id: {}".format(self.snapId, firstHaloId))
                 
            print(nextProg, lastsnap.groupLocs) 
            if nextProg in lastsnap.groupLocs:
                print("a;sdlkfja;eijoeiwafhasdf")
                yLoc = DC(redshifts[lastsnap.groupLocs[nextProg][1]])*sToYr/c
                plt.plot([self.xs[ii], lastsnap.groupLocs[nextProg][0]], [DC(redshifts[self.snapId])*sToYr/c,yLoc], 'k')
            ii += 1



    def remove_small_groups(self, limit):
        removalList = []
        for group in self.fofGroups:
            if self.fofGroups[group].m < limit:
                removalList.append(group)

        # remove all of the necessary groups from the list
        for kk in removalList:
            del self.fofGroups[kk]





