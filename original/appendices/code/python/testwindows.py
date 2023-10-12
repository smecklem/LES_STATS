
from numpy import linspace,diff,mean,arange,set_printoptions,zeros
from pylab import plot,show

def getwindows(times, windowtime, nwindows):
    n = times.size
    maxtime = times[-1]
    
    windowlen = int(n*windowtime/maxtime)
    centerdiff = n//(nwindows+1)
    ws = windowlen//2
    if windowlen>n: raise Exception("windowtime too high!")
    
    #print('maxtime:', maxtime)
    #print('nwindows: ', nwindows)
    #print('windowtime: ',windowtime)
    #print('windowlen: ',windowlen)
    #print('centerdiff: ',centerdiff)
    #print('ws: ',ws)
    
    center = [(i+1)*centerdiff for i in range(nwindows)]
    
    # suggested starts and ends
    ss = [max(i-ws, 0)  for i in center]
    es= [min(i+ws+1,n) for i in center]

    starts = []
    ends = []
    lens = []
    for si,ei in zip(ss,es):
        sl= si
        el= si+windowlen

        sh = ei - windowlen
        eh = ei

        if el<=n:
            starts.append(sl)
            ends.append(el)
            lens.append(el-sl)
        else:
            starts.append(sh)
            ends.append(eh)
            lens.append(eh-sh)

    assert all([i==lens[0] for i in lens])
    return starts,ends

def plotwindows(times, v, start, end):
    up = mean(v)/5.0
    for s,e in zip(start, end):
        plot(times[s:e], v[s:e] + up)
        up = up*-1.0
    show()
    return

if __name__=="__main__":
    set_printoptions(linewidth = 200)
    n = 20
    v =  arange(10.0,n+10.0)
    times = arange(0.0, n)
    nwindows = 3
    windowtime = 5.5
    dt = mean(diff(times))

    start, end = getwindows(times, windowtime, nwindows)
    for s,e in zip(start, end): print(s,e,e-s)

    testarray = zeros((nwindows+2,n))
    testarray[0,:] = arange(n)
    testarray[1,:] = v
    for i,s,e in zip(range(nwindows), start, end):
        testarray[i+2,s:e] = 1
    print(testarray)

    plotwindows(times, v, start, end)
