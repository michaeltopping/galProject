import numpy as np






suffix = "br.fits"

filename = "../shapley2003_spec/all.cat"
mappingfile = "../spec/mask_design/z.ssa22a.current"
outfile = "../shapley2003_spec/allnew.cat"
out = open(outfile, 'w')
data = np.genfromtxt(filename, dtype=str, filling_values="xxx", usecols=(0,1,2,3,4))
mapdata = np.genfromtxt(mappingfile, dtype=str)
for row in data:
    filename = row[1]
    name = row[0]
    zem = row[2]
    zabs = row[3]
    note = row[4]

    if name in mapdata[:,2]:
        index = np.where(name==mapdata[:,2])[0][0]
        zem = mapdata[index,8]
        zabs = mapdata[index,9]
    out.write("{:13s}    {}    {}    {}    {}\n".format(name, filename, zem, zabs, note))    
print(mapdata[:,1])
