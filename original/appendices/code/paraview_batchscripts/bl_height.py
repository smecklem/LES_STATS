"""
Batch compute boundary layer heights? 

@author: Nick Gibbons
"""

from paraview.simple import *
from glob import glob
from numpy import zeros,radians,array,cos,sin,tan,arctan,linspace,nanmax,argmax,unique,argmin,isnan,save
from pylab import plot,show,title,ylabel,xlabel
from scipy.interpolate import interp1d

filename = 'stats.vtk'
soln = OpenDataFile(filename)

def get_array_from_points(points, arrayname):
    adata = points.GetArray(arrayname)
    n = adata.GetDataSize()
    a = zeros(n)
    adata.ExportToVoidPointer(a)
    return a
    
# Cone parameters from Tristan's CAD
theta = arctan((68.0-62.0)/(130.0-70.0)) # Cone Half Angle
C = array([750.0, 0.0, 0.0])/1000.0      # Cone point location
ahat = array([-1.0, 0.0, 0.0])           # Cone axis normal vector
omega = radians(15.0)

# Generate startting and ending points in the injector symmetry plane coords
xstartcoord = 0.0925
xendcoord = 0.142
N = 20
L = 10.0e-3
if N==3: raise Exception("Possible broadcast error inc, proceed with caution")

# Vectors in the injector symmetry plane reference frame
xdashed = linspace(-1.0*(C[0]-xstartcoord), -1.0*(C[0]-xendcoord), N)
ydashed = xdashed*tan(theta) # Careful with quadrants here.
zdashed = 0.0*xdashed

spd = zeros((3,N))
spd[0,:] = xdashed
spd[1,:] = ydashed
spd[2,:] = zdashed
wnd = array([-sin(theta), cos(theta), 0.0]) # Wall normal vector
wpd = array([cos(theta), sin(theta), 0.0]) # Wall perpendicular vector, aligned with flow toward C

print("spd", spd[:,0])
print("wnd ",wnd [:])
print("wpd",wpd[:])

# convert to rotated frame
romega = array([[1.0, 0.0, 0.0],[0.0, cos(omega),-sin(omega)],[0.0,sin(omega),cos(omega)]])
sp = romega.dot(spd)
wn = romega.dot(wnd)
wp = romega.dot(wpd)

sp = (sp.T).copy()
ep = sp+L*wn
print("sp", sp[0])
print("wn", wn)

# translate to actual frame
s = sp+C
e = ep+C
print("start[0]", s[0])
print("end[0]", e[0])
print("Loading Data")
Show(soln)
temps = []
vels = []
resolution = 200

for i in range(N):
    print("Getting data for line",i)
    srt = s[i]
    end = e[i]
    line = PlotOverLine(Input=soln,Source='High Resolution Line Source')
    line.Source.Point1 = srt 
    line.Source.Point2 = end
    line.Source.Resolution = resolution
    Show(line)

    # Paraview's arc-length is apparently bugged so i'll do it myself smh
    #dvec = end - srt
    #linepoints = srt + outer(dvec,linspace(0,1,resolution)).T
    #arc = linspace(0,L,resolution+1) # Whyyyyyyyyy + 1?

    data = servermanager.Fetch(line)
    points = data.GetPointData()

    u = get_array_from_points(points, 'um')
    v = get_array_from_points(points, 'vm')
    w = get_array_from_points(points, 'wm')
    vel = u*wp[0] + v*wp[1] + w*wp[2]
    vels.append(vel)

    Delete(line)
Delete(soln)

vels = array(vels)
arc = linspace(0,L,resolution+1) 
blhs = []
for x,vel in zip(s[:,0],vels):
    start = argmin(isnan(vel))
    end = vel.size - argmin(isnan(vel[::-1]))

    veltrimmed = vel[start:end]
    htrimmed = arc[start:end]

    velunique,idxs = unique(veltrimmed, return_index=True)
    hunique = htrimmed[idxs]

    velint = interp1d(hunique,velunique)
    n = 2000
    nh = linspace(hunique[0], hunique[-1], n)
    nvel = velint(nh)

    #maxvel = nanmax(nvel)
    maxvel = 2200.0
    isbl = nvel>0.99*maxvel
    blidx = argmax(isbl)
    if blidx==0: blidx = nh.size-1
    blh = nh[blidx] - nh[0]

    print("nunique: ", velunique.size)
    print("blh: ", blh)
    #plot(nh, nvel, 'b.')
    plot(nh, nvel, 'b-')
    plot(hunique, velunique,'r.')
    plot(nh[blidx], nvel[blidx],'go')
    title("Velocity profile @ x={}".format(x))
    show()

    blhs.append(blh)

blhs = 1000*array(blhs)
blxs = 1000*s[:,0]
save('scripts/npydata/blxs',blxs)
save('scripts/npydata/blhs',blhs)
save('scripts/npydata/vels',vels)

title("Boundary Layer Height: golf (99% {} m/s)".format(maxvel))
ylabel('BL Height (mm)')
xlabel("x position (mm)")
plot(blxs, blhs)
show()

# Note that there may be some nans in the lineplots becase the lines could be starting or
# ending out of the domain
#print(s[12], e[12])
#for i,v in enumerate(vels):
#    #maxvel = nanmax(v) # FIXME: Unsafe! Line might go out of domain before reaching the edge of the BL
#    maxvel = 2200.0 # yeah you need a better method of detecting the boundary layer end.
#    isbl = v>0.99*maxvel # May be some warnings here, could conceal problems
#    blidx = argmax(isbl)
#    blh = arc[blidx]
#    blv = v[blidx]
#    blhs.append(blh)
#    if 0<i and i<8:
#        plot(v,arc)
#        plot(blv, blh, 'ro')
#show()
#
#print(len(blhs), sp.shape)
#plot(s[:,0], blhs)
#show()

print "... Done"

