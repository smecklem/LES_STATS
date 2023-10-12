"""
Convergence plot for hotel

@author: Nick Gibbons 
"""

from libmonitor import parse
from pylab import plot,show,semilogx,loglog,semilogy,legend,close,savefig,title,xlabel,ylabel
from numpy import array, correlate, mean, sum, diff, fft, zeros
from testwindows import getwindows, plotwindows
from glob import glob

tf = 2.80700977895e-05 # from compute_flowtime.py

files = sorted(glob('../monitors/monitor.*.dat'))
#files = ['../logs/monitor.05.dat']

for file in files:
    data = parse(file)
    print(file, data.keys())

    offset_times = data['time'] - data['time'][0]
    startup_time = 3.0*tf
    latetimes = offset_times>startup_time

    times = offset_times[latetimes]
    n = times.size
    dt = mean(diff(times))
    varname = 'p'

    var = data[varname][latetimes]
    varf = var - mean(var)
    show()

    # Split into windows
    nwindows = 3
    windowtime = 3*tf
    start,end = getwindows(times, windowtime, nwindows)
    #plotwindows(times, var, start, end)

    nsignal = end[0]-start[0] # All windows should be the same length
    freq = fft.rfftfreq(nsignal, d=dt)
    print("freq.shape", freq.shape)
    meanspec = zeros(nsignal//2+1)

    for s,e in zip(start, end):
        spec = fft.rfft(varf[s:e])
        print("spec.shape", spec.shape)
        meanspec+=abs(spec)

    meanspec/=float(nwindows)

    number = file.split('.')[-2]
    print("THING", number)
    loglog(freq, meanspec, label=varname)
    xlabel('Hz')
    title('Pressure Fluctuation Spectrum: {}'.format(number))
    savefig('../vis/spectra/{}.svg'.format(varname+'-'+number))
    close()
    show()

print("Done")
#show()


