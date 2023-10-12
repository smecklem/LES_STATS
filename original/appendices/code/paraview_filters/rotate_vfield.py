
"""
Compute rotated velocity field

@author: Nick
"""
from paraview.numpy_support import vtk_to_numpy
import vtk.numpy_interface.dataset_adapter as dsa
import vtk.numpy_interface.algorithms as algs
import numpy as np
from numpy import zeros, einsum, minimum


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

u = get_array('u')
v = get_array('v')
w = get_array('w')

theta = radians(-15.0) # First set of injectors
udash = u
vdash = v*cos(theta) - w*sin(theta)
wdash = v*sin(theta) + w*cos(theta)

post_array(udash, 'udash')
post_array(vdash, 'vdash')
post_array(wdash, 'wdash')

