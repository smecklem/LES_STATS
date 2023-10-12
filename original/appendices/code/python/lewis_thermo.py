"""
Module to lookup species properties in the NASA Lewis thermodynamic database.

@author: Nick Gibbons
"""

from string import ascii_letters
from collections import Counter
from numpy import log,linspace,array
from numpy import sum as npsum

letters = set(ascii_letters)
DBPATH = '/home/qungibbo/programs/ceq/thermo.inp'

def readdb(name):
    """ Retrieve species 'name' from the .db file, using a naive approach """
    with open(DBPATH) as fp:

        iterline=iter(fp) # Make an iterator for inner loop behaviour

        # Start looking for the beginning of a species entry
        for line in iterline:
            if not line[0] in letters: continue # Skip if not a letter in first position

            if line.startswith(name):        # This could be the one, check the rest
                
                header = line.strip().split()
                if header[0]!=name: continue # Nope false alarm

                lines = [line]               # We've found it!
                for nextline in iterline:
                    if nextline[0] in letters: break
                    lines.append(nextline)
                break                        # Break the outer for loop and jump to return lines
        else:
            raise Exception("Name: {} not found!".format(name))
    return lines

def get_species(name):
    """ Pointless rename to avoid clashes with 'species' variables """
    return species(name)

class species(object):
    """ Container for species data """
    def __init__(self, name):
        data = readdb(name)
        header = data[0].strip().split()
        if header[0]!=name:
            raise IOError("Database read failed! {}!={}".format(name, header[0]))

        self.name = name
        self.ref = ' '.join(header[1:])
        info = data[1].strip().split()
        self.intervals = int(info[0])
        
        # FixMe: Do the atomic constiuents thing by column
        self.Hf298 = float(info[-1])     # J/mol
        self.M     = float(info[-2])/1e3 # kg/mol
        self.hf298 = self.Hf298/self.M      # J/kg

        self.ai = []
        self.bi = []
        self.Ti = []
        self.ji = []
        for i in range(self.intervals):
            hdr= data[2+3*i]
            Tlow, Thigh = hdr[:22].strip().split()
            therest = hdr[24:].strip().split()
            self.Ti.append((float(Tlow),float(Thigh))) # Temperature bounds
            self.ji.append([float(j) for j in therest[:7]])
            self.dH298 = float(therest[-1]) # Weird place for H0(298.0) - H0(0)

            a = data[3+3*i].replace('D','E') # Fix Fortran Double nonsense
            afloat = [float(a[16*j:16*(j+1)]) for j in range(5)]

            a = data[4+3*i].replace('D','E') # Fix Fortran Double nonsense
            afloat.extend([float(a[16*j:16*(j+1)]) for j in range(2)])
            self.ai.append(afloat) 

            self.bi.append([float(a[16*j:16*(j+1)]) for j in range(3,5)])

        self.ai = array(self.ai)
        self.bi = array(self.bi)
        self.Ti = array(self.Ti)
        self.ji = array(self.ji)

        self.Hf = self.Hf298 #FIXME: This is wrong before, it's an approximation now
        self.hf = self.Hf/self.M
        self.R = 8.3144598/self.M # (J/mol/K) Universal Gas Constnat
        
        self.cp0 = self.cp(298.15)
        self.atoms = self.countatoms(self.name)
        self.natoms = sum(self.atoms.values())
        self.cv0 = self.cp0 - self.R

        # QM degrees of freedom, divided by two. 
        self.cf = 1.5
        if self.natoms==1:
            self.cf += 0.0
        elif self.natoms==2:
            self.cf += 1.0
        elif self.natoms>2:
            self.cf = self.cp(200.1)/self.R - 1.0
        return

    #def Trange_index(self,T):
    #    """ Return the index needed for a specific temperature lookup """
    #    print(self.Ti)
    #    for i,Tvals in enumerate(self.Ti):
    #        Tlow, Thigh = Tvals
    #        if Tlow<T and T<=Thigh:
    #            return i
    #    raise Exception("No Data for Temperature: {} ({}-{})".format(T, self.Ti[0][0], self.Ti[-1][-1]))

    def Trange_index(self,T):
        """ Return the index needed for a specific temperature lookup """
        lb = []
        ub = []
        and_ = []
        #xor_ = []
        
        for Tvals in self.Ti:
            Tlow, Thigh = Tvals
            l = Tlow<T
            u = T<=Thigh
            lb.append(l)
            ub.append(u)
            and_.append((1*l)*(1*u))
            #xor_.append(((1*l)+(1*u))%2)
        
        check = 1-sum(i for i in and_) # Check for entries with all bounds missing
        overflow = check*(1*lb[0])*(len(lb)-1) # Low ones get zero, add len-1 to get high ones
        index = sum(i*j for i,j in enumerate(and_))
        index = index+overflow
        return index

    def cp(self, T):
        """ Return the empircal specific heat capacity @ const. pressure """
        idx = self.Trange_index(T)
        #return sum(a*T**j for a,j in zip(self.ai[idx], self.ji[idx]))*self.R # J/kg/K
        #cp = npsum(self.ai[idx]*(T**self.ji[idx].T).T,axis=-1)*self.R
        Cp0onR = self.Cp0onR(T)
        cp = Cp0onR*self.R # Note that the Cp0onR returns molar Cp/Ru, convert here to cp (J/kg)
        return cp

    def Cp0onR(self,T):
        """ Return Cp0/R, equation 4.9 in nasacea_I """
        #Cp0onRT = sum(a*T**j for a,j in zip(self.ai[idx], self.ji[idx])) 
        idx = self.Trange_index(T)
        a0,a1,a2,a3,a4,a5,a6 = self.ai[idx].T
        a7,a8 = self.bi[idx].T

        Cp0onR  =  a0*T**-2
        Cp0onR +=  a1/T
        Cp0onR +=  a2
        Cp0onR +=  a3*T
        Cp0onR +=  a4*T**2
        Cp0onR +=  a5*T**3
        Cp0onR +=  a6*T**4
        return Cp0onR

    def H0onRT(self,T):
        """ Return H0/RT, equation 4.10 in nasacea_I """
        idx = self.Trange_index(T)
        a0,a1,a2,a3,a4,a5,a6 = self.ai[idx].T
        a7,a8 = self.bi[idx].T

        H0onRT = -a0*T**-2
        H0onRT+=  a1*log(T)/T
        H0onRT+=  a2
        H0onRT+=  a3*T/2.0
        H0onRT+=  a4*T**2/3.0
        H0onRT+=  a5*T**3/4.0
        H0onRT+=  a6*T**4/5.0
        H0onRT+=  a7/T
        return H0onRT

    def S0onR(self,T):
        """ Return S0/RT, equation 4.11 in nasacea_I """
        idx = self.Trange_index(T)
        a0,a1,a2,a3,a4,a5,a6 = self.ai[idx].T
        a7,a8 = self.bi[idx].T

        S0onR  = -a0*T**-2/2.0
        S0onR += -a1/T
        S0onR +=  a2*log(T)
        S0onR +=  a3*T
        S0onR +=  a4*T**2/2.0
        S0onR +=  a5*T**3/3.0
        S0onR +=  a6*T**4/4.0
        S0onR +=  a8
        return S0onR

    def countatoms(self,name):
        """ Set up the atomic constituents count for this species """
        numbers = set('123456789')
        lowers = set('abcdefghijklmnopqrstuvwxyz')
        uppers = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ') 

        # First split into Actual Elements (HHe2-> H, He2)
        element = name[0] # Assume first character is a capital
        elements = []
        for char in name[1:]:
            if char not in uppers:
                element+=char
            else:
                elements.append(element)
                element = char

        elements.append(element)

        atoms = []
        # Now deal with numbers
        for element in elements:
            atom = ''.join([char for char in element if char not in numbers])

            if len(atom)==len(element)>0:
                num = 1
            else:
                num = int(element[len(atom):]) # Extract the number from the end
            atoms.extend([atom]*num)

        return Counter(atoms)


    def __repr__(self):
        return '\n'.join(map(str,list(self.__dict__.items())))

                
if __name__=='__main__':
    # Basic Test:
    H = species('H')
    print('H: \n', H)
    H2 = species('H2')
    print('H2: \n', H2)
    H2O= species('H2O')
    print('H2O: \n',H2O)
    print("H2.cp(T)",H2.cp(300))
    print("H2.cp(linspace(300,4000,3))",H2.cp(linspace(300,4000,3)))
    print("H2.H0onRT(T)",H2.H0onRT(300))
    print("H2.H0onRT(linspace(300,4000,3))",H2.H0onRT(linspace(300,4000,3)))
    
