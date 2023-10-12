#!/usr/bin/python3
"""
Monitor files from US3D are formatted so they can be plotted in Tecplot
but since I'm a dinosaur I'm going to use Python instead. This is exactly
the same as the postmonitor file in the /utils directory but this one
is accesable from PYTHONPATH

TODO: the executable postmonitor script should just call this one :thinking:?

@author: Nick Gibbons
"""

#from pylab import plot,show,xlabel,ylabel,legend,title
from matplotlib import pyplot as plt
import matplotlib
from glob import glob
from copy import copy
from argparse import ArgumentParser
from string import ascii_letters as letters
from itertools import cycle
from numpy import array, bool, zeros, ndarray

COLOURS = ['b','g','r','c','m','y','k']
LINES = ['-','--','-.',':']
MARKERS = ['o', 'v', '*', 'h']
matplotlib.rcParams['svg.fonttype'] = 'none'
LETTERS = set(letters)

def cell_pos(file):
    """ Retreive the cell center coordinates from a US3D monitor file"""
    with open(file) as fp:
        for line in fp:
            if line.startswith('#   Nearest match ='):
                foundline = line.strip()
                break
        else:
            raise Exception("No nearest match line!")

    x = float(foundline[-12*3:-12*2].replace('D','e'))
    y = float(foundline[-12*2:-12*1].replace('D','e'))
    z = float(foundline[-12*1:].replace('D','e'))
    return array([x,y,z])


def parse(filename,start=0):
    """ Extract a Dictionary of things from a monitor file """
    with open(filename) as fp:
        lines = [line for line in fp]

    iterlines = iter(lines)
    name = filename.replace(".dat","")

    # Due to annoying line breaks, we need to be clever about
    # deciding what is header and what isn't
    header = ''
    for line in iterlines:
        if line.startswith('#'): 
                continue
        # Quick count of all nonblank number entries
        linechars = line.strip().replace(' ','')
        lineletters= sum(1 for i in linechars if i in LETTERS)
        if lineletters>len(linechars)/2:
            # This line has mostly letters and we're probably still in header
            header += line.strip()
        else: 
            # We're through the header, so branch and make a new list
            numbers = [line.strip().split()]
            break

    # Annoying header splitting problem fix
    for line in iterlines:
        numbers.append(line.strip().split())

    # Now load the numbers into a list of lists
    numbers = zip(*numbers) # Transpose rows to columns
    numbers = [array(list(map(makenumber,i))) for i in numbers]
    numbers = [d[start:].copy() for d in numbers]
    for d in numbers: print(d.shape)


    # Finally parse the reassembled header
    header = header.split('=')[1].strip().split()

    if len(numbers)!=len(header):
        print(" #### WARNING #### ")
        print("Data Length Mismatch in file: ", filename)
        print(" ")

    data = dict(zip(header, numbers))
    data['name'] = name
    return data

def fix_restart_glitches(data,verbose=True):
    """ Fix reversals in time from restarts in the CFD """
    iters = data['iter']
    print("Starting with {} entries".format(iters.size))

    keep = (zeros(iters.shape)+1).astype(bool)
    lastiter = iters[0] - 1
    scrubbed = 0

    for i in range(iters.size):
        iter = iters[i]

        if iter!=lastiter+1:
            if verbose: print("Found jump @ iter", iter)
            for j in range(1,i+1):
                previter = iters[i-j]
                keep[i-j] = False
                scrubbed += 1
                if verbose: print("Scrubbed iter: ", previter)
                if previter==iter:
                    lastiter = iter
                    if verbose: print("Finished scrubbing\n", lastiter)
                    break
        else:
            lastiter = iter

    print("Found {} datapoints to remove".format(scrubbed))
    newdata = {}
    for k,v in data.items():
        if type(v)==ndarray:
            newdata[k] = v[keep].copy()
        else:
            newdata[k] = v

    newsize = newdata['iter'].size
    print("New size is {}, diff is {}".format(newsize, iters.size-newsize))
    return newdata


def makefloat(instring):
    """ Handle some glitches in the floating points """
    try:
        return float(instring)
    except(ValueError):
        return 0.0


def makenumber(instring):
    """ Handle some glitches in the incoming numbers """
    try:
        return int(instring)
    except(ValueError):
        try:
            return float(instring)
        except(ValueError):
            return 0.0
    return 0.0
        

def makeplot(filedata, plotnames, plottitle=''):
    """ Plot each variable in plotnames, with a line from each
        case in filedata.
    """
    print("Begin Make Plot: ")
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colours = cycle(iter(COLOURS))

    for d in filedata:
        markers = cycle(iter(MARKERS))
        lines   = cycle(iter(LINES))
        colour = next(colours)

        for plotvar in plotnames:
            if not plotvar in d.keys():
                print("No data for {} in {}".format(plotvar,d['name']))
                continue
            if len(d[plotvar])>14:
                line   = next(lines)
                fmt = colour+line 
            else:
                marker = next(markers)
                fmt = colour+marker          
            labelstring = d['name']+":"+plotvar+""
            print("Making plot {} with fmt: {}".format(labelstring, fmt))
            ax.plot(d['time'], d[plotvar],fmt,label=labelstring)

    ax.set_xlabel('Time: sec')
    ax.set_ylabel(plotvar)
    fig.suptitle(plottitle)
    ax.tick_params(axis='x')
    ax.tick_params(axis='y')
    ax.legend(loc=0)
    plt.show()
    return
    
def init():
    """ Collect files matching the input arguments """
    print("US3D Monitor File Plotting tool: ")
    parser = ArgumentParser(description="Plot Monitor Files")
    parser.add_argument("-f", "--files", type=str, nargs='+',
                      help="Space delimited filenames")

    parser.add_argument("-v", "--variables", type=str, nargs='+',
                      help="Space delimited variables to plot")
    parser.add_argument("-n", "--normalize", action='store_true',
                        help="Whether to normalize data")
    parser.add_argument("-t", "--title", type=str, nargs=1,
                        help="Plot title")
    args = parser.parse_args()

    variables = args.variables
    expfilenames = []
    for name in args.files:
        expfilenames.extend(glob(name))
   
    if len(expfilenames)==0:
        raise Exception("No file(s) found that match input args")

    if args.title==None: args.title = ['']
    
    print("Files Found: ")
    for i in expfilenames: print(i)
    print(" ")

    return expfilenames, variables, args

def normalize(filedata):
    """ Normalize each entry by its first entry. """
    for data in filedata:
        for k in data.keys():
            if k=='name': continue
            if k=='time': continue
            firstentry = copy(data[k][0])
            data[k] = [t/firstentry for t in data[k]]
    return filedata

if __name__=='__main__':
    files, variables, args = init()
    filedata = [parse(name) for name in files]
    if args.normalize: normalize(filedata)
    makeplot(filedata,variables,args.title[0]) 






