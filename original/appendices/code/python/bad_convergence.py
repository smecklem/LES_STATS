"""
Zach has suggested I add create some "bad" plots to illustrate the concept.

@author: Nick Gibbons
"""

from glob import glob
from libmonitor import *
import matplotlib.pyplot as plt
from numpy import linspace, cumsum, arange, zeros, mean, diff, floor, ceil, max, sin, min, cos
from scipy.interpolate import interp1d

plt.rcParams['svg.fonttype'] = 'none'
import matplotlib as mpl

if __name__=='__main__':
    ft = 2.80700977895e-05 # Assume the same as hotel 
    #files = sorted(glob('../oldlogs/monitor.*.dat'))
    files = ['../monitors/monitor.01.dat',
             '../monitors/monitor.02.dat',
#             '../monitors/monitor.03.dat',
             '../monitors/monitor.04.dat',
             '../monitors/monitor.05.dat']
    datas = [parse(file) for file in files]
    datas = [fix_restart_glitches(data, verbose=False) for data in datas]
    for d in datas: d['flowtimes'] = (d['time'] - d['time'][0])/ft

    fig = plt.figure(figsize=(12,10))
    axess = [fig.add_subplot(2,2,i+1) for i in range(len(datas))]

    #plotthing = 'u'
    plotthing = 'p'

    if plotthing=='u':
        vkey  = 'u'
        vlabel= 'u'
        vunits = 'm/s'
    elif plotthing=='p':
        vkey  = 'p'
        vlabel= 'p'
        vunits = 'Pa'
    else:
        raise Exception("Wrong plotthing! ({})".format(plotthing))

    tt = datas[2]['flowtimes']
    pt = datas[2][plotthing]
    nt = tt.size

    step = zeros(pt.shape)
    step[nt//2:] += 0.3*mean(pt)
    ttbad = linspace(8,10,tt.size)
    ptbad = interp1d(tt,pt)(ttbad)

    badsignals = [
        (tt, ptbad),
        (tt,sin(tt*47)*(max(pt) - min(pt)) + mean(pt) + 0.5*pt),
        (tt,pt + sin(tt/8)*0.5*(max(pt) - min(pt))),
        (tt,pt + step)
    ]

    badnames = [
        'No Short Wavelengths       ',
        'Overly Regular       ',
        'Long wavelength oscillations       ',
        'Step Change       '
    ]

    for file,data,axes in zip(badnames,badsignals,axess):
        times, thing = data
        axes.plot(times, thing, label = file)

        axes.set_ylabel('{} ({})'.format(vkey, vunits))
        axes.set_xlabel('(t-t0)/tf')
        axes.legend(loc=0)
    plt.tight_layout()
    plt.show()

    print("...Done")
