"""
Compute the change in wall temperature using a semi-infinite transient analysis

@author: Nick
"""

from paraview.numpy_support import vtk_to_numpy
import vtk.numpy_interface.dataset_adapter as dsa
import vtk.numpy_interface.algorithms as algs
import numpy as np
from numpy import zeros


print("Begin...")
def get_array(name):
    for inp in inputs:
        if name in inp.CellData.keys():
            arr = inp.CellData[name]
            return arr
    else:
        raise KeyError('{} not found'.format(name))

def numpy2vtk(arr,dset,aa):
    vtkdata = dsa.numpyTovtkDataArray(arr)
    vtkarray = dsa.vtkDataArrayToVTKArray(vtkdata,dset)
    vtkarray.Association = aa
    return vtkarray

def post_array(arr,name):
    print("Saving", name, "...")
    output.CellData.append(arr, name)
    return

randomarray = inputs[0].CellData[inputs[0].CellData.keys()[0]] # Yes, this happens.
numpy_to_vtk = lambda arr: numpy2vtk(arr, randomarray.DataSet, randomarray.Association)

# get paraview.vtk.dataset_adapter.VTKArray object
print "Getting arrays to numpy"

# Aluminium properties (SI)
k = 205.0    # W/m2
rho = 2710.0 # kg/m3
c = 910.0    # J/kg/K
q = get_array('qw')
texposure = 5e-3 # approximate test article expose time from pitot survey shot 12262, condition A

dT = 2*q*np.sqrt(texposure)/np.sqrt(np.pi*rho*c*k)
post_array(dT,'dT')

print "... Done"
