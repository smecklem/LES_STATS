"""
Test code for extracting cells in a slice file that are along a boundary:

Notes; I can't beleive this actually worked...
@author: Nick
"""

import h5py
from numpy import arange, array, where, in1d, save

# Target face 5
fstring = "5      3      5  faces  246797297  247283696     486400        *wall"

zone,itype,zoneid,ztype,startf,endf,nf,name = fstring.split()

startf = int(startf)
endf   = int(endf)
nf     = int(nf)

print(zone,itype,zoneid,ztype,startf,endf,nf,name)

faces = arange(startf-1, endf)
print(faces.size, nf, faces.size==nf)

print("Reading files...")
try:
    slices = h5py.File('../slices.h5','r')
    slicecells = slices['definitions/cells/list_1'][:]
    print(slicecells[0:10], slicecells.size)
finally:
    slices.close()

try:
    conn = h5py.File('../conn.h5','r')
    ife = conn['ife'][:]
    facecells= ife[faces,0]
    print(facecells[0:10], facecells.size)
finally:
    conn.close()

# Strangely this actually seems to be the best way to do this:
print("Computing set intersection...")
scset = set(slicecells)
fcset = set(facecells)
intersection = scset.intersection(fcset)
sfcells = array(sorted(list(intersection)))
print(sfcells[0:10], sfcells.size)
save('npydata/dcells_face{}.npy'.format(zoneid), sfcells)

# So sfcells are the global indexes of the cells in the slice dataset that are also on the target face
# You need to convert these global indexes to their positions in the slice array
sliceindexes = where(in1d(slicecells,sfcells))[0]
print(sliceindexes[0:10], sliceindexes.size)

print("Saving...")
save('npydata/si_face{}.npy'.format(zoneid, sliceindexes), sliceindexes)
#save('npydata/ransi_face{}.npy'.format(zoneid, sliceindexes), sliceindexes)
print("...Done")
