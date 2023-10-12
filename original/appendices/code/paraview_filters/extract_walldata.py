"""
Compute wall heat transfer from slice data.

@author: Nick
"""

from paraview.numpy_support import vtk_to_numpy
import vtk.numpy_interface.dataset_adapter as dsa
import vtk.numpy_interface.algorithms as algs
import numpy as np
from numpy import zeros, einsum, minimum, array


def get_array(name):
    for inp in inputs:
        if name in inp.CellData.keys():
            arr = inp.CellData[name]
            return arr
    else:
        raise KeyError('{} not found'.format(name))

def post_array(arr,name):
    print("Saving", name, "...")
    output.CellData.append(arr, name)
    return

randomarray = inputs[0].CellData[inputs[0].CellData.keys()[0]] # Yes, this happens.

nx = randomarray.size
print("Found arrays, size: {}".format(nx))

try:
    print("Loading grid co-oords's @ cell centers")
    xce= load('scripts/npydata/slice_cell_centers.npy')
except(IOError):
    print("Computing grid co-oords's @ cell centers")
    xce= zeros((nx,3))
    for i in range(nx):
        xci = zeros(3)
        vtk_cells = inputs[0].GetCell(i) 
        npoints = vtk_cells.GetNumberOfPoints()
        points = vtk_cells.GetPoints()
        for j in range(npoints):
            point = points.GetPoint(j)
            xci += array(point)
        xci/=npoints
        xce[i,:] = xci.copy()
    save('scripts/npydata/slice_cell_centers.npy', xce)

# Cone parameters from Tristan's CAD
theta = np.arctan((68.0-62.0)/(130.0-70.0)) # Cone Half Angle
C = array([750.0, 0.0, 0.0])/1000.0      # Cone point location
ahat = array([-1.0, 0.0, 0.0])           # Cone axis normal vector

# See log entry for 23/07/18
A = xce - C                    
Amag = np.sqrt(A[:,0]**2 + A[:,1]**2 + A[:,2]**2)
phi = np.arccos(A.dot(ahat)/Amag)  
psi = theta - phi 
pmag = Amag*np.sin(psi)
Emag = Amag*np.cos(psi)

faceindexes = load('scripts/npydata/si_face5.npy')
iswall = zeros(nx)
iswall[faceindexes] = 1.0

Tf = get_array('t')[faceindexes]
ktf= get_array('kt')[faceindexes]
pf = get_array('p')[faceindexes]
Emagf = Emag[faceindexes]
pmagf = pmag[faceindexes]
xf = xce[faceindexes,0]
yf = xce[faceindexes,1]

dx = 1e-4/1000.0/2.0 # Half of the first cell height, from hotel_grid/cluspec.txt, converted from mm
Tw = 300.0
wht = ktf*(Tf - Tw)/dx

save('scripts/npydata/xf.XXX.npy',   xf)
save('scripts/npydata/wht.XXX.npy',  wht)
save('scripts/npydata/pf.XXX.npy', pf)
save('scripts/npydata/pmagf.XXX.npy', pmagf)

#save('scripts/npydata/RANS.xf.npy',   xf)
#save('scripts/npydata/RANS.wht.npy',  wht)
#save('scripts/npydata/RANS.pf.npy', pf)
#save('scripts/npydata/RANS.pmagf.npy', pmagf)

#post_array(iswall, 'iswall')
#post_array(pmag, 'pmag')
#post_array(Amag, 'Amag')
#post_array(Emag, 'Emag')

