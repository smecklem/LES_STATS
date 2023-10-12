"""
Plot data to show statistical convergence of simulation "golf"

@author: Nick Gibbons
"""

from glob import glob
from libmonitor import *
import matplotlib.pyplot as plt
from numpy import linspace, cumsum, arange, zeros, mean, diff, floor, ceil
from scipy.interpolate import interp1d

plt.rcParams['svg.fonttype'] = 'none'
import matplotlib as mpl

def cume_av_plot(fig, axes, data, vkey, vlabel):
    thing = data[vkey]
    cav = cumsum(thing)/(arange(thing.size) + 1)
    ncrs = data['iter']

    axes.plot(ncrs, thing, label='${}$'.format(vlabel))
    axes.plot(ncrs, cav, 'k-', linewidth=2)

    return

def reverse_av_plot(fig, axes, data, vkey, vlabel):
    print("final time: ", data['flowtimes'][-1])
    times = data['flowtimes']
    thing = data[vkey]
    ncrs = data['iter']
    timestep = data['dt']

    pmm = (cumsum(thing[::-1]*timestep[::-1])/cumsum(timestep[::-1]))[::-1]

    axes.plot(times, thing, label='${}$'.format(vlabel))
    axes.plot(times, pmm, 'k-', linewidth=2)

    return

def add_stats_lines(fig, axes, data, stats_iters):
    ymin, ymax = axes.get_ylim()
    ncrs = data['iter']
    flowtimes = data['flowtimes']

    interpolator = interp1d(ncrs,flowtimes,fill_value=(flowtimes.min(), flowtimes.max()))
    stats_times = [interpolator(it) for it in stats_iters]

    for time in stats_times:
        axes.plot([time, time], [ymin, ymax],'k--',linewidth=2)
    return

def add_iters(fig, axes, data):
    raise Exception("This is broken: FIXME!")
    time = data['time']
    ncrs = data['iter']
    flowtimes = data['flowtimes']

    timeticks = axes.get_xticks()
    print(timeticks)

    interpolater = interp1d(flowtimes,ncrs,fill_value=(ncrs.min(), ncrs.max()),bounds_error=False)
    iterticks = interpolater(timeticks)
    print(iterticks)

    axes2= axes.twiny()
    axes2.set_xticklabels([str(int(i)) for i in iterticks])
    axes2.set_xlabel(r'iter')
    return axes2

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

    stats_iters = [6001,16001,32000]
    for file,data,axes in zip(files,datas,axess):
        idx = file.split('.')[-2]
        reverse_av_plot(fig, axes, data, vkey, vlabel+'_{r}'.replace('r',idx))

        add_stats_lines(fig, axes, data, stats_iters)

        #add_iters(fig, axes, data)
        axes.set_ylabel('{} ({})'.format(vkey, vunits))
        axes.set_xlabel('(t-t0)/tf')
        axes.legend(loc=0)
    plt.tight_layout()
    plt.show()

    print("...Done")
