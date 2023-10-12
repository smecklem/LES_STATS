"""
Script File for paraview console

@author: Nick Gibbons
"""

from paraview.simple import *
from glob import glob
from numpy import sin,cos,radians

print "Begin..."
OFFSET = 0.0001 #So we aren't inside the domain
angle = radians(5.09)
nxo = 0.0
nyo = -sin(angle)
nzo = cos(angle) 
vxo = OFFSET*nxo
vyo = OFFSET*nyo
vzo = OFFSET*nzo
names = glob('monitors/monitor*.dat')
names.sort()
assert len(names)>0


for filename in names:
    with open(filename) as fp:
       for line in fp:
           if 'match' in line:
               numbers = line.strip().split(' = ')[1]
               numbers = numbers.replace('D','e')
               break
       else:
           raise Exception('No Line Found')
    l = len(numbers)/3
    x = float(numbers[0*l:1*l])
    y = float(numbers[1*l:2*l])
    z = float(numbers[2*l:3*l])

    pointname = '.'+filename.split('.')[1]
    text = a3DText(Text = pointname)
    SetDisplayProperties(text, Scale=[0.001, 0.001, 0.001])
    SetDisplayProperties(text, Origin=[x+vxo,y+vyo,z+vzo])
    SetDisplayProperties(text, Orientation=[16.0, 0.0, 0.0])
    SetDisplayProperties(text, DiffuseColor=[1.0, 1.0, 1.0])
    
    Show(text)
Render()
print "...Done"
