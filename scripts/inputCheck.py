import numpy as np






suffix = "br.fits"

filename = "../shapley2006_spec/all.cat"
mappingfile = "./ssa22a.oldmaskkey.txt"
outfile = "../shapley2006_spec/allnew.cat"
out = open(outfile, 'w')
data = np.genfromtxt(filename, dtype=str, filling_values="xxx", usecols=(0,1,2,3,4))
mapdata = np.genfromtxt(mappingfile, dtype=str)
for row in data:
    filename = row[0]+"."+row[1]+"."+suffix
    name = row[1]
    zem = row[2]
    zabs = row[3]
    note = row[4]

    if filename[0:-5] in mapdata[:,1]:
        index = np.where(filename[0:-5]==mapdata[:,1])[0][0]
        name = mapdata[index,0][7:]
        zem = mapdata[index,5]
        zabs = mapdata[index,6]
        note = mapdata[index,2]
    out.write("{:13s}    {}    {}    {}    {}\n".format(name, filename, zem, zabs, note))    
print(mapdata[:,1])
