#### import the simple module from the paraview
from paraview.simple import *
from numpy import sin,cos,radians
#### disable automatic camera reset on 'Show'
# find source
solvtk = FindSource('sol.vtk')

# create a new 'Slice'
slice1 = Slice(Input=solvtk)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# Properties modified on slice1.SliceType
L = 0.03
x0 = 0.1
#theta = radians(12.501) # First set of injectors
theta = radians(5.0) # First set of injectors
#theta = radians(17.5) # Central symmetry plane
#theta = radians(20.0) # second set of injectors
#theta = radians(22.499) # Far wall
origin = [x0,-L*cos(theta),-L*sin(theta)]
print(origin)
normal = [0.0,-sin(theta), cos(theta)]
print(normal)

slice1.SliceType.Origin = origin
slice1.SliceType.Normal = normal

#slice1.SliceType.Origin = [0.08, -0.0289777747887, -0.00776257135308]
#slice1.SliceType.Normal = [0.0, -0.258819045102982, 0.965925826288945]

Hide(solvtk)
Show(slice1)
