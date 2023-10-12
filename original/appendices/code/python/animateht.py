"""
How to animate arrays in a portable way? Not Like This that's for sure.

@author: Nick Gibbons
"""

from glob import glob
from numpy import load,argsort,zeros,linspace,array,mean,sort,diff,argmax
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from intervalmaker import get_90percent_CI
from us3dutils import scrape_attrs, humantime
import matplotlib as mpl
mpl.style.use('classic')

startedtime = 6.944905794E-03 # from plot_convergence.py

xfs= [load(file) for file in sorted(glob('npydata/xf.*.npy'))]
whts= [load(file) for file in sorted(glob('npydata/wht.*.npy'))]
rawtimes = scrape_attrs('../slices.h5', 'results/cells/list_1', 'time')
rawtimes = array(list(map(float,rawtimes)))
rawtimes = sort(rawtimes)

inc = rawtimes>startedtime
print("Rawtimes: {} Include only: {}",format(len(rawtimes), str(sum(inc))))

times = rawtimes - rawtimes[0]
htimes = array([humantime(t) for t in times])
dt = mean(diff(times))
dt = humantime(dt)
print("Mean dt", dt)
print("len(whts)", len(whts))

idxs = argsort(xfs[0])
xf = [i[idxs] for i in xfs]
wht = [i[idxs] for i in whts]

assert len(xf)==len(wht)
assert len(xf)==len(htimes)
assert len(xf)>0
assert len(wht)>0

xf = array(xf)
wht = array(wht)
mwht = mean(wht, axis=0)
lwht, uwht = get_90percent_CI(wht)
print("lwht.shape", lwht.shape, type(lwht))

# Normalise everything for Tristan
mpiastart = 0.1
xf -= mpiastart

m_to_mm = 1000.0
xf *= m_to_mm 

wht_unit = 100000.0
wht /= wht_unit
mwht /= wht_unit
uwht /= wht_unit
lwht /= wht_unit

images = []
fig = plt.figure(figsize=(10,6))
axes = plt.gca()
axes.set_xlim([xf[0][0], xf[0][-1]])
axes.set_ylim([0.0, 5])
#axes.set_ylim([min(wht[0]), max(wht[0])])
#axes.set_aspect(70.0/0.4*500.0/1046.34) # To match the visualisations in mixing_0*.svg
plt.xlabel('x-from-mpia (mm)')
plt.ylabel('WHT (UNITS)')
plt.title("Wall Heat Transfer Data (dt={})".format(dt))
plt.tight_layout()

print(len(xf), len(htimes))
for i in range(len(xf)):
    thisimage = []
    thisimage.extend(axes.plot(xf[i],wht[i],'g-',linewidth=1))
    thisimage.extend(axes.plot(xf[0],mwht,  'g--',linewidth=1))
    thisimage.append(axes.annotate('t={}'.format(htimes[i]),xy=(10, 10), xycoords='figure pixels'))
    images.append(thisimage)

fill = plt.fill_between(xf[0],lwht, uwht, alpha=0.2, facecolor='green')
plt.legend([images[0][0], images[0][1]], labels=['LES HT', 'Mean HT'])

ani = animation.ArtistAnimation(fig, images, interval=200, repeat_delay=0,blit=True) 
#plt.show()
ani.save('test.mp4',bitrate=2**11)
