"""
How to animate arrays in a portable way? Not Like This that's for sure.

@author: Nick Gibbons
"""

from glob import glob
from numpy import load,argsort,zeros,linspace,array,mean,sort,diff,sum
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from intervalmaker import get_90percent_CI
from us3dutils import scrape_attrs, humantime
import matplotlib as mpl
mpl.style.use('classic')


def get_wht(xfs, whts, rawtimes, startedtime):
    """ Compute mean and variances for plotting wht """
    inc = rawtimes>startedtime
    print("Rawtimes: {} Include only: {}".format(len(rawtimes), str(sum(inc))))
    rawtimes = rawtimes[inc]
    xfs = [xfs[i]  for i,istrue in enumerate(inc) if istrue]
    whts= [whts[i] for i,istrue in enumerate(inc) if istrue]

    times = rawtimes - rawtimes[0]
    htimes = array([humantime(t) for t in times])
    dt = mean(diff(times))
    dt = humantime(dt)
    print("Mean dt", dt)
    print("len(whts)", len(whts), len(xfs),len(htimes))

    idxs = argsort(xfs[0])
    xf = [i[idxs] for i in xfs]
    wht = [i[idxs] for i in whts]

    assert len(xf)==len(wht)
    assert len(xf)==len(htimes)
    assert len(xf)>0
    assert len(wht)>0

    xf = array(xf)[0].copy()
    wht = array(wht)
    mwht = mean(wht, axis=0)
    lwht, uwht = get_90percent_CI(wht)
    print("lwht.shape", lwht.shape, type(lwht))
    d = {'xf':xf, 'mwht':mwht, 'lwht':lwht, 'uwht':uwht}
    return d

# hotel numbers
startedtime = 0.01191099603784225 # from plot_convergence
xfs = [load(file) for file in sorted(glob('/media/qungibbo/data/sims/mpia/hotel/scripts/npydata/xf.*.npy'))]
whts= [load(file) for file in sorted(glob('/media/qungibbo/data/sims/mpia/hotel/scripts/npydata/wht.*.npy'))]
rawtimes = scrape_attrs('/media/qungibbo/data/sims/mpia/hotel/slices.h5', 'results/cells/list_1', 'time')
rawtimes = array(list(map(float,rawtimes)))
rawtimes = sort(rawtimes)
dhotel = get_wht(xfs, whts, rawtimes, startedtime)

# golf numbers
startedtime = 6.944905794E-03 # from plot_convergence.py
xfs = [load(file) for file in sorted(glob('./npydata/xf.*.npy'))]
whts= [load(file) for file in sorted(glob('./npydata/wht.*.npy'))]
print(len(xfs))
rawtimes = scrape_attrs('../slices.h5', 'results/cells/list_1', 'time')
rawtimes = array(list(map(float,rawtimes)))
rawtimes = sort(rawtimes)
dgolf = get_wht(xfs, whts, rawtimes, startedtime)

images = []
fig = plt.figure(figsize=(10,6))
axes = plt.gca()
axes.set_xlim([dgolf['xf'][0], dgolf['xf'][-1]])

axes.set_ylim([0.0, 5e5])
plt.xlabel('x (m)')
plt.ylabel('WHT (W/m2)')
plt.title("Wall Heat Transfer Data Comparison")

thisimage = []
thisimage.extend(axes.plot(dgolf['xf'],dgolf['mwht'], 'g-',linewidth=1))
thisimage.extend(axes.plot(dhotel['xf'],dhotel['mwht'], 'r-',linewidth=1))
images.append(thisimage)

fill = plt.fill_between(dgolf['xf'],dgolf['lwht'], dgolf['uwht'], alpha=0.2, facecolor='green')
fill = plt.fill_between(dhotel['xf'],dhotel['lwht'], dhotel['uwht'], alpha=0.2, facecolor='red')
plt.legend([images[0][0], images[0][1]], labels=['Double MPIA HT','Single MPIA HT'])

plt.tight_layout()
plt.show()
