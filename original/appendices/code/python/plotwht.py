"""
How to animate arrays in a portable way? Not Like This that's for sure.

@author: Nick Gibbons
"""

from glob import glob
from numpy import load,argsort,zeros,linspace,array,mean,sort,diff,sum,argmax
import matplotlib.pyplot as plt
from intervalmaker import get_90percent_CI
from us3dutils import scrape_attrs, humantime, readcsv
#plt.switch_backend('agg')
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

# golf numbers
startedtime = 6.944905794E-03 # from plot_convergence.py
xfs = [load(file) for file in sorted(glob('./npydata/xf.*.npy'))]
whts= [load(file) for file in sorted(glob('./npydata/wht.*.npy'))]
print(len(xfs))
rawtimes = scrape_attrs('../slices.h5', 'results/cells/list_1', 'time')
rawtimes = array(list(map(float,rawtimes)))
rawtimes = sort(rawtimes)
dgolf = get_wht(xfs, whts, rawtimes, startedtime)

# rans golf numbers
data = readcsv('npydata/ransline.csv')
xfrans = data['Points:0']
whtrans = data['qw']
idxs = argsort(xfrans)
xfrans = xfrans[idxs]
whtrans = whtrans[idxs]

# Adjust longfoxtrot data to use foxtrot conditions
xflam = load('/media/qungibbo/data/sims/mpia/longfoxtrot/scripts/npydata/xf.RANS.npy')
whtlam = load('/media/qungibbo/data/sims/mpia/longfoxtrot/scripts/npydata/wht.RANS.npy')
idxs = argsort(xflam)
xflam = xflam[idxs]
whtlam = whtlam[idxs]

xflam2 = load('/media/qungibbo/data/sims/mpia/foxtrot/scripts/npydata/xf.XXX.npy')
whtlam2 = load('/media/qungibbo/data/sims/mpia/foxtrot/scripts/npydata/wht.XXX.npy')
idxs = argsort(xflam2)
xflam2 = xflam2[idxs]
whtlam2 = whtlam2[idxs]

data = readcsv('/media/qungibbo/data/sims/mpia/longfoxtrot/scripts/npydata/turb_slice.csv')
xfturb = data['Points:0']
whtturb = data['qw']
idxs = argsort(xfturb)
xfturb = xfturb[idxs]
whtturb = whtturb[idxs]

i2 = argmax(xflam2>dgolf['xf'][0])
i = argmax(xflam>dgolf['xf'][0])
diffht = whtlam[i] - whtlam2[i2]
whtlam -= diffht
whtturb-= diffht
print("sizes", whtlam.size, whtturb.size)

# Normalise everything for Tristan
mpiastart = 0.1
dgolf['xf'] -= mpiastart
xfrans -= mpiastart 
xflam -= mpiastart

m_to_mm = 1000.0
dgolf['xf'] *= m_to_mm 
xfrans *= m_to_mm 
xflam *= m_to_mm 

wht_unit = 100000.0
dgolf['mwht'] /= wht_unit
dgolf['lwht'] /= wht_unit
dgolf['uwht'] /= wht_unit
whtturb /= wht_unit
whtrans /= wht_unit
whtlam /= wht_unit

images = []
fig = plt.figure(figsize=(10,6))
axes = plt.gca()
axes.set_xlim([dgolf['xf'][0], dgolf['xf'][-1]])

axes.set_ylim([0.0, 5])
#plt.xlabel('x (m)')
#plt.ylabel('WHT (W/m2)')
plt.xlabel('x-from-MPIA (mm)')
plt.ylabel('WHT (UNITS)')
plt.title("Wall Heat Transfer Data Comparison")


thisimage = []
thisimage.extend(axes.plot(dgolf['xf'],dgolf['mwht'], 'g-',linewidth=1))
thisimage.extend(axes.plot(xfrans,whtrans, 'g--',linewidth=1))
#thisimage.extend(axes.plot(xflam,whtlam, 'k--',linewidth=1))
#thisimage.extend(axes.plot(xfturb,whtturb, 'k.',linewidth=1))
images.append(thisimage)

fill2= plt.fill_between(xflam,whtlam, whtturb[:-3], alpha=0.2, facecolor='black')
fill1= plt.fill_between(dgolf['xf'],dgolf['lwht'], dgolf['uwht'], alpha=0.2, facecolor='green')
plt.legend([images[0][0], images[0][1], fill2], labels=['LES HT','RANS HT','No Injection HT'],loc='upper left')

#plt.legend([images[0][0], images[0][1], images[0][2]], labels=['LES HT','RANS HT','Laminar HT'],loc='upper left')


plt.tight_layout()
#plt.savefig('asdf.png')
#plt.close()
plt.show()
