"""
Batch compute boundary layer heights? 

@author: Nick Gibbons
"""

from paraview.simple import *
from glob import glob
from numpy import zeros,radians,array,cos,sin,tan,arctan,linspace,nanmax,argmax,unique,argmin,isnan,save
from pylab import plot,show,title,ylabel,xlabel
from scipy.interpolate import interp1d

def get_array_from_points(points, arrayname):
    adata = points.GetArray(arrayname)
    n = adata.GetDataSize()
    a = zeros(n)
    adata.ExportToVoidPointer(a)
    return a
    
def setup_sweep_points(xstart, xend, L, N=20):
    # Cone parameters from Tristan's CAD
    theta = arctan((68.0-62.0)/(130.0-70.0)) # Cone Half Angle
    C = array([750.0, 0.0, 0.0])/1000.0      # Cone point location
    ahat = array([-1.0, 0.0, 0.0])           # Cone axis normal vector
    omega = radians(15.0)

    # Generate startting and ending points in the injector symmetry plane coords
    if N==3: raise Exception("Possible broadcast error inc, proceed with caution")

    # Vectors in the injector symmetry plane reference frame
    xdashed = linspace(-1.0*(C[0]-xstart), -1.0*(C[0]-xend), N)
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
    return wp,s,e

def get_vel_from_line(soln, srt, end, wp, resolution=200):
    line = PlotOverLine(Input=soln,Source='High Resolution Line Source')
    line.Source.Point1 = srt 
    line.Source.Point2 = end
    line.Source.Resolution = resolution
    Show(line)

    data = servermanager.Fetch(line)
    points = data.GetPointData()

    u = get_array_from_points(points, 'um')
    v = get_array_from_points(points, 'vm')
    w = get_array_from_points(points, 'wm')
    vel = u*wp[0] + v*wp[1] + w*wp[2]
    Delete(line)
    return vel

def get_high_resolution_bl_profile(arc,vel,n=2000):
    start = argmin(isnan(vel))
    end = vel.size - argmin(isnan(vel[::-1]))
    if end!=vel.size:
        print("WARNING: DOMAIN CLIPPING!!!!!!!!!!")
        print("Start", start, "end", end)

    veltrimmed = vel[start:end]
    htrimmed = arc[start:end]

    velunique,idxs = unique(veltrimmed, return_index=True)
    hunique = htrimmed[idxs]

    velint = interp1d(hunique,velunique)
    nh = linspace(hunique[0], hunique[-1], n)
    nvel = velint(nh)
    return nh, nvel
    
def get_bl_height(nh,nvel,criterion):
    isbl = nvel>criterion
    blidx = argmax(isbl)
    if blidx==0: blidx = nh.size-1
    blh = nh[blidx] - nh[0]
    blv = nvel[blidx]
    return blh,blv

print("Loading Data")
filename = 'stats.vtk'
soln = OpenDataFile(filename)
Show(soln)

N = 40
resolution = 200
xstart = 0.0935
xend = 0.142
L = 10.0e-3
n = 2000
arc = linspace(0,L,resolution+1) 
wp,s,e = setup_sweep_points(xstart, xend, L, N)

vels = []
blhs = []
blxs = []
for i in range(N):
    print("Getting data for line",i)
    srt = s[i]
    end = e[i]
    vel = get_vel_from_line(soln, srt, end, wp) 
    vels.append(vel)

    nh,nvel = get_high_resolution_bl_profile(arc,vel,n)

    #maxvel = nanmax(nvel)
    maxvel = 2200.0
    criterion = 0.99*maxvel
    blh,blv = get_bl_height(nh,nvel,criterion)
    blhs.append(blh)
    blxs.append(srt[0])

    print("blh", blh)
    if i==4:
        plot(nh, nvel, 'b-')
        plot(blh+nh[0], blv,'go')
        title("Velocity profile @ x={}".format(srt[0]))
        show()
        vels.pop()
        blhs.pop()
        blxs.pop()

blhs = 1000*array(blhs)
blxs = 1000*array(blxs)
blus = array(vels)
save('scripts/npydata/blxs',blxs)
save('scripts/npydata/blhs',blhs)
save('scripts/npydata/blus',vels)

plot(blxs, blhs, 'b-')
xlabel('x (mm)')
ylabel('h (mm)')
title("Boundary Layer Heights u={}".format(criterion))
show()

Delete(soln)

print "... Done"

