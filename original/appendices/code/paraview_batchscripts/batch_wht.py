"""
Batch compute wht data

@author: Nick Gibbons
"""

from paraview.simple import *
from glob import glob

filenames = sorted(glob('./slices/slices.s0001*.vtk'))
print("Found nfiles: ", len(filenames))

with open('./scripts/extract_walldata.py') as fp:
     script = fp.read()

for i,name in enumerate(filenames):
    print "Computing file: ", name, '{}/{}'.format(i+1,len(filenames))
    soln = OpenDataFile(name)
    prgfil = ProgrammableFilter(soln)
    

    my_script = script.replace('XXX',str(i).zfill(3))
    prgfil.Script = my_script
    Show(prgfil)
    Delete(prgfil)
    Delete(soln)

print "... Done"

