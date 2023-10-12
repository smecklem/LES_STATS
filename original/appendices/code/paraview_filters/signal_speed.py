"""
Compute timestep requirements for LES

@author: Nick
"""
import sys
print("exe: ", sys.executable)

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
volume = algs.volume(inputs[0])
L = np.cbrt(volume)

a = get_array('a')
u = get_array('u')
v = get_array('v')
w = get_array('w')

vel = (u**2 + v**2 + w**2)**0.5
signalspeed = vel+a
signaltime = L/signalspeed

post_array(volume, 'volume')
post_array(signalspeed, 'signalspeed')
post_array(signaltime*1e9, 'signaltime(ns)')
post_array(L, 'L')

