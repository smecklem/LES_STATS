"""
Compute viscosity of a gas misture using the Blottner curvefits

Based on US3D's gastrans

@author: Nick Gibbons
"""

from numpy import array,zeros,sqrt,log,exp
from lewis_thermo import get_species

DBPATH = '/home/qungibbo/programs/us3d/1.0-RC22.12/props/blottner.db'

def make_blott_from_species_list(species):
    blottner_db = readdb(DBPATH)

    blot = zeros((3,len(species)))
    for i,name in enumerate(species):
        blot[:,i] = blottner_db[name]
    return blot

def readdb(name):
    """ Read all the species in file 'name' """
    data = {}
    with open(name) as fp:

        iterline=iter(fp) # Make an iterator for inner loop behaviour
        for line in iterline:
            if line.strip().startswith('[BLOTTNER]'):
                break

        for line in iterline:
            if line.strip().startswith('[/BLOTTNER]'):
                break
            if line.strip().startswith('!'):
                continue
            line = line.replace(',','')
            name,A,B,C = line.strip().split()
            name = name.strip("'") 
            A = float(A.replace('d','e'))
            B = float(B.replace('d','e'))
            C = float(C.replace('d','e'))
            data[name] = [A,B,C]

    if len(data)==0: raise Exception("No items found!")
    return data 

def gastrans(blot,M,rs,tt):
    """ Compute the viscosity of a gas mixture """
    ns = len(rs)

    rmrs = zeros((ns,ns))
    pin  = zeros((ns,ns))
    for n in range(ns):
       for m in range(ns):
          pin[n,m]  = 1.0/sqrt(8.0*(1.0 + M[n]/M[m]))
          rmrs[n,m] = (M[n]/M[m])**2.50e-1

    t1 = log(tt)
    rms = zeros(rs.shape)
    chi = zeros(rs.shape)
    for n in range(ns):
       rms[n] = 0.1*exp((blot[0,n]*t1 + blot[1,n])*t1 + blot[2,n])
       chi[n] = rs[n]/M[n]

    rmi = 1.0/rms
    rmm = chi.sum(axis=0)
    chi = chi/rmm

    rmt = 0.0
    for n in range(ns):
       phi = 0.0
       for m in range(ns):
          t1  = 1.0 + sqrt(rms[n]*rmi[m])*rmrs[m,n]
          phi = phi + chi[m]*t1*t1*pin[n,m]
       t1 = rms[n]*chi[n]/phi      
       rmt = rmt + t1
    rmu = rmt
    return rmu

def gastrans2(blot,gassp,rs,tt):
    """ Compute the viscosity of a gas mixture """
    M = array([i.M for i in gassp])
    mu = gastrans(blot,M,rs,tt)
    return mu

def gastrans2_kt(blot,gassp,rs,tt,vib=False):
    """ Compute the thermal conductivity of a gas mixture """
    Ru = 8.3144621 # J/K/mol (REEEEEEEEEEEEEEEEEEEE)
    M = array([i.M for i in gassp])
    cf = array([i.cf for i in gassp])
    cvv = array([i.cp(tt) - i.cp0 for i in gassp])

    ns = len(rs)
    Rs = Ru/M

    rmrs = zeros((ns,ns))
    pin  = zeros((ns,ns))
    for n in range(ns):
       for m in range(ns):
          pin[n,m]  = 1.0/sqrt(8.0*(1.0 + M[n]/M[m]))
          rmrs[n,m] = (M[n]/M[m])**2.50e-1

    t1 = log(tt)
    rms = zeros(rs.shape)
    chi = zeros(rs.shape)
    for n in range(ns):
       rms[n] = 0.1*exp((blot[0,n]*t1 + blot[1,n])*t1 + blot[2,n])
       chi[n] = rs[n]/M[n]

    rmi = 1.0/rms
    rmm = sum(chi)
    chi = chi/rmm

    rkt = 0.0
    rkv = 0.0
    for n in range(ns):
       phi = 0.0
       for m in range(ns):
          t1  = 1.0 + sqrt(rms[n]*rmi[m])*rmrs[m,n]
          phi = phi + chi[m]*t1*t1*pin[n,m]
       t1 = rms[n]*chi[n]/phi      
       rkt = rkt + t1*(cf[n] + 2.25)*Rs[n]
       rkv = rkv + t1*cvv[n]

    if vib: rkt = rkt + rkv
    return rkt

def test_mu():
    blottner_db = readdb('/media/qungibbo/data/sims/dns/bravo/h2-air13_jach92.dat')
    for i in blottner_db.items(): print(i)
    spnames = 'N2,O2,H2,H2O,OH,HO2,H2O2,NO,NO2,HNO,N,H,O'.split(',')
    test_header = ["vtkOriginalPointIds","p","rho","t","tv","m","csN2","csO2","csH2","csH2O","csOH","csHO2","csH2O2","csNO","csNO2","csHNO","csN","csH","csO","u","v","w","a","mu","vtkOriginalCellIds","vtkValidPointMask","arc_length","Points:0","Points:1","Points:2"]
    test_line = [164117,1.1075e+05,0.023241,1641.1,0,0.23051,0.24267,4.6151e-05,0.67094,0.08314,3.8044e-05,1.1602e-08,1.6545e-10,3.0299e-07,1.6893e-15,4.8942e-09,3.4345e-10,0.0031613,3.3646e-07,490.75,-277.84,191.2,2583.2,3.4417e-05,33651,1,0,0.013579,0.0099666,-0.0022819]

    test = dict(zip(test_header,test_line))
    ns = len(spnames)
    species = [get_species(i) for i in spnames]
    M = array([i.M for i in species])

    testrs = array([test['cs'+i] for i in spnames])*test['rho']

    blot= zeros((3,ns))
    for i,sp in enumerate(spnames):
        blot[:,i] = blottner_db[sp]

    mu = gastrans(blot, M, testrs, test['t'])
    mu2= gastrans2(blot,species, testrs, test['t'])
    print("mu", mu)
    print("mu2", mu2)
    print('testmu', test['mu'])
    print(" ")

def test_kt():
    blottner_db = readdb('/media/qungibbo/data/sims/mpia/india/gas.dat')
    for i in blottner_db.items(): print(i)
    spnames = 'N2,O2,C2H4'.split(',')
    test_header = ["vtkOriginalPointIds","rho","u","v","w","t","p","m","csN2","csO2","csC2H4","sa-mut","kt","mu","vtkOriginalCellIds","vtkValidPointMask","arc_length","Points:0","Points:1","Points:2"]
    test_line = [279154,0.016372,1220.2,140.4,38.97,923.03,4361.4,2.1176,0.66066,0.22022,0.11912,0.00046535,0.06397,3.7508e-05,130066,1,3.6102e-05,0.14909,-0.056845,-0.015252]


    test = dict(zip(test_header,test_line))
    ns = len(spnames)
    species = [get_species(i) for i in spnames]
    M = array([i.M for i in species])

    testrs = array([test['cs'+i] for i in spnames])*test['rho']

    blot= zeros((3,ns))
    for i,sp in enumerate(spnames):
        blot[:,i] = blottner_db[sp]

    mu = gastrans(blot, M, testrs, test['t'])
    mu2= gastrans2(blot,species, testrs, test['t'])
    kt = gastrans2_kt(blot,species,testrs,test['t'],vib=True)
    print("mu", mu)
    print("mu2", mu2)
    print('testmu', test['mu'])
    print(" ")

    print('kt', kt)
    print('testkt', test['kt'])
    print(" ")

def test_kt_2D():
    blottner_db = readdb('/media/qungibbo/data/sims/mpia/india/gas.dat')
    for i in blottner_db.items(): print(i)
    spnames = 'N2,O2,C2H4'.split(',')
    test_header = ["vtkOriginalPointIds","rho","u","v","w","t","p","m","csN2","csO2","csC2H4","sa-mut","kt","mu","vtkOriginalCellIds","vtkValidPointMask","arc_length","Points:0","Points:1","Points:2"]
    test_line1= [279154,0.016372,1220.2,140.4,38.97,923.03,4361.4,2.1176,0.66066,0.22022,0.11912,0.00046535,0.06397,3.7508e-05,130066,1,3.6102e-05,0.14909,-0.056845,-0.015252]
    test_line2= [245547,0.045229,2187.6,268.05,72.08,336.95,4382.4,5.9934,0.75,0.25,9.5175e-67,1.2749e-07,0.0288,2.1065e-05,96492,1,0.00025703,0.12996,-0.041515,-0.011223]
    test_line = array([test_line1, test_line2]).T

    test = dict(zip(test_header,test_line))
    print(test)
    ns = len(spnames)
    species = [get_species(i) for i in spnames]

    testrs = array([test['cs'+i] for i in spnames])*test['rho']

    blot= zeros((3,ns))
    for i,sp in enumerate(spnames):
        blot[:,i] = blottner_db[sp]

    mu2= gastrans2(blot,species, testrs, test['t'])
    kt = gastrans2_kt(blot,species,testrs,test['t'],vib=True)
    print("mu2", mu2)
    print('testmu', test['mu'])
    print(" ")

    print('kt', kt)
    print('testkt', test['kt'])
    print(" ")

if __name__=='__main__':
    #test_mu()
    #test_kt()
    test_kt_2D()
