"""
Check logs

@author: Nick Gibbons
"""

from glob import glob
from us3dutils import extract_datafiles, stitch_datafiles
from pylab import plot,show,ylabel,xlabel,semilogy,title
import matplotlib as mpl
mpl.style.use('classic')

def print_dict_to_table(d):
    size = max(len(k) for k in d.keys())
    ksorted = sorted(d.keys())
    for k in ksorted:
        print(k.rjust(size), ":", d[k])

files = glob('../logs/us3d*.log')

ft = 2.80700977895e-05
mds = extract_datafiles(*files)

print("")
print("Found {} files\n".format(len(mds)))

for i,file in enumerate(mds):
    print(i,'.)')
    d = {}

    for key in file.keys():
        if key!='data':
            d[key] = file[key]
    try:
        d['flowtimes'] = (d['final_time'] - d['start_time'])/ft
    except(KeyError):
        pass 
    print_dict_to_table(d)
    print(" ")

rans = [mds[3]]
data = stitch_datafiles(rans)

title("RANS Stage Convergence")
ylabel(r"$ \sqrt{\Delta \rho \Delta \rho_i}$")
xlabel(r"Iteration Number")
semilogy(data['iter'], data['residual'])
show()
