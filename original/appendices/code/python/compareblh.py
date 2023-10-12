"""
Compare longfoxtrot and golf boundary layer stuff.

@author: Nick Gibbons
"""

from numpy import load
from pylab import plot,show,ylabel,xlabel,title,legend,ylim

#gx = load('npydata/blxs.npy')
#gh = load('npydata/blhs.npy')
#
#fx = load('/media/qungibbo/data/sims/mpia/foxtrot/scripts/npydata/blxs.npy')
#fh = load('/media/qungibbo/data/sims/mpia/foxtrot/scripts/npydata/blhs.npy')

lfx =load('/media/qungibbo/data/sims/mpia/longfoxtrot/scripts/npydata/blxs.npy')
lfh =load('/media/qungibbo/data/sims/mpia/longfoxtrot/scripts/npydata/blhs.npy')

ix = load('/media/qungibbo/data/sims/mpia/india/scripts/npydata/blxs.npy')
ih = load('/media/qungibbo/data/sims/mpia/india/scripts/npydata/blhs.npy')


mpiastart = 0.095202736777608*1000
mpiaend = 0.104861437900961*1000

title("Boundary Layer Thickness: MPIA Experimental Model")
ylabel('BL Height (mm)')
xlabel('Distance from Leading Edge (mm)')
#plot(gx, gh, label='LES + MPIA')
plot(ix, ih, label='LES + MPIA')
plot(lfx, lfh, label='Laminar No MPIA')
plot([mpiastart,mpiastart], [0,8.5],'k--')
plot([mpiaend,mpiaend], [0,8.5],'k--')
legend(loc='lower left')
ylim(0,8.5)
show()
