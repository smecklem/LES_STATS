"""
(Allegedly) Organised Library of helpful stuff for working with US3D.

@author: Nick Gibbons
"""

import h5py
import numpy as np
import argparse
from glob import glob
from copy import copy

def scrape_attrs(filename, path, attr):
    """ This function scrapes every object in the group specified by 'path'
    and makes a list of any attributes they have that are called 'attr'. Note  
    this method requires some knowledge of how the file is structured: If you
    don't know or care, try the old version in ~/Sim../laser../test_in../res../
    or '$ h5ls -v -r myfile.h5' to get information about the file structure.

    >> scrape_attrs('filename.h5', 'results/cells/list_1', 'time')

    @author: Nick
    """
    fp = h5py.File(filename, 'r')
    outlist = []
    try:
        group = fp[path]

        for thing in list(group.values()):
            for name, value in thing.attrs.items():
                if name==attr: 
                    if len(value)==1: outlist.append(value[0])
                    else: outlist.append(value)

    finally:
        fp.close()

    return outlist

def us3d_file(filename):
    """ Extract detailed solution information from filename. """
    fp = h5py.File(filename,'r')
    #print("Details of {}/info/solver".format(fp.filename))
    things = {}
    try:
        group = fp['info/solver']
        for k,v in group.attrs.items():
            things[k] = v
    finally:
        fp.close()

    return things

def stitch_datafiles(mdlist):
    """ Take the output of extract data and generate a single 'data' dictionary that
    is just the columns from each logfile, merged into a single chronological order
    with any unsaved sections removed.

    In: md list -- Expects list(dict.keys() - > [name, stop_iter, etc... , data] )
    """
    data = {}

    for md in mdlist:
        print("Extending with:  {} {} to: \n".format(md['name'], md['start_iter']))
        print("                 {} {}".format(' '*len(md['name']), md['stop_iter']))
        for k,v in md['data'].items():
            slicefromthisfile = v[:md['stop_iter']-md['start_iter']+1]
            data.setdefault(k,list()).extend(slicefromthisfile)
    return data

def extract_datafiles(*filepatterns):
    """
    Python routine for extracting US3D residual data from log files of stdout.
    You should only use this if you don't want to use convergence.dat 
    for some reason. The actual logic is handled by parse_lines, which may be
    a little memory intensive for large files. This is the main control 
    function, which controls agglomerating data for multiple files, specified 
    using Unix-like wilcarding. Example usage:

    data = extract_data('us3dJul07-*.log') # Pool data from July 7th

    @author: Nick
    """
    filenames = []
    for filepattern in filepatterns:
        filenames.extend(glob(filepattern))

    if len(filenames)==0:
        raise Exception("No Match For {}".format(filepatterns[0]))

    metadata = []
    filedata = []
    for filename in filenames:
        print("Reading file: {}...".format(filename))
        with open(filename) as fp:
            lines = [line.strip() for line in fp]
        print("    -- {} lines".format(len(lines)))
        if len(lines)<=1:
            print("No Lines... Skipping")
            continue

        md = parse_metalines(lines)
        if 'date' not in list(md.keys()):
            print("    -- Read Error! (No date)")
            continue
        if 'time' not in list(md.keys()):
            print("    -- Read warning! (No time)")
        md['name'] = filename
        metadata.append(md)
        filedata.append(parse_lines(lines))

    # Now sort both lists by the date in the metadata one
    zipdata = list(zip(metadata, filedata))
    zipdata.sort(key = lambda x: x[0]['date'])
    metadata, filedata = list(zip(*zipdata))

    for md, fd in zip(metadata, filedata):
        md['data'] = fd # Add the actual data to the md dict here

    return list(metadata) # Not just meta data anymore

def parse_lines(lines):
    """
    Parse a list of strings corresponding to lines in a file opened
    somewhere else. We want to extract the convergence info from a
    us3d runlog, keeping in mind that there may be varying amounts
    of header and footer stuff, as well as breaks in the numbers
    where stdout logs writing etc.

    @author: Nick
    """
    # We need two for loops. The prettiest way to do this uses a named
    # iterator, which is looped over directly.
    lineiterator = iter(lines) 
    header = ''
    for line in lineiterator:
        if line.startswith('iter '):
            header = line.split()
            next(lineiterator) # Skip the --- -----   ---- stuff
            break
    else:
        return dict(zip(list(), list())) # No iter lines found...
        
    
    # We're now through the header, record anything that looks
    # like convergence info, by checking iteration number.
    try: numbers = [next(lineiterator).split()]
    except(StopIteration): 
        return dict(zip(header, []*len(header)))

    iterno = int(numbers[0][0])

    for line in lineiterator:
        if line.startswith(str(iterno+1)+" "):
            numbers.append(line.split())
            iterno+=1

    # Convert strings to numbers (int or float) and tranpose with zip(*)
    numbers = [numcast(list(i)) for i in zip(*numbers)] 

    data = dict(zip(header,numbers))
    return data



def numcast(inlist):
    """ Convert a list of strings to numbers (if possible) based on 
    the inferred type of the first element. Items that fail int and
    float conversions should be returned as they were. 

    @author: Nick
    """
    value = inlist[0]
    castf = lambda x: x
    try:
        int(value)
        castf = int
    except ValueError:
        try:
            float(value)
            castf = float
        except ValueError:
            castf = lambda x: x
    return list(map(castf,inlist))
       
def extract_metadata(filenames):
    """ We would like to agglomerate some metadata from logfiles. This
    includes number of cores, number of iterations, and walltime.

    @author: Nick
    """
    #filenames = []
    #for filepattern in filepatterns:
    #    filenames.extend(glob(filepattern))

    if len(filenames)==0:
        raise Exception("No Matches For {}".format(str(filepatterns)))

    data = []
    aborted = 0

    for filename in filenames:
        with open(filename) as fp:
            lines = [line.strip() for line in fp]
        datum = parse_metalines(lines)
        datum['name'] = filename

        if 'time' not in list(datum.keys()):
            aborted+=1 
            print("No Time found for file: {}".format(filename))
            continue
        if 'date' not in list(datum.keys()):
            aborted+=1
            continue
        data.append(datum)

    if aborted == 1:
        print((" Warning: {} log indicates an aborted simulation.".format(aborted)))
    if aborted > 1:
        print((" Warning: {} logs indicate aborted simulations.".format(aborted)))
    return data

def parse_metalines(lines):
    """ Parse an individual files metadata """
    datum = {}
    lineiterator = iter(lines)
    unitconv={'secs':1,'mins':60,'hrs':3600}
    monconv=[None,'Jan','Feb','Mar','Apr','May','Jun',
             'Jul','Aug','Sep','Oct','Nov','Dec']
    teqtable = {'SA-Catris':1, 'SA-neg':1, 'DES97':1,'DDES':1, 'IDDES':1, 'SST-V':2}
    datum['neq'] = 4 # default for ru,rv,rw,rE
    datum['chemistry'] = False
    datum['turbulence'] = None
    datum['machine'] = 'unspecified'

    for line in lineiterator:
        if 'Current date: ' in line:
            strings = line.strip().split()
            month = strings[4]
            month = str(monconv.index(month)).zfill(2)
            day = strings[5].zfill(2)
            time = strings[6]
            year = strings[7]
            dstamp = '.'.join([year,month,day]) + '-' + time
            datum['date'] = dstamp
        if 'a part of US3D' in line:
            version = line.strip().split()[4]
            datum['version'] = version
        if line.startswith('-- Number of processors'):
            datum['np'] = int(line.split()[-1])
        if '-- Number of cells:' in line:
            ncells= line.strip().split()[4]
            ncells = ncells.replace(',','')
            datum['ncells'] = int(ncells)
        if '-- Number of faces:' in line:
            nfaces= line.strip().split()[4]
            nfaces = nfaces.replace(',','')
            datum['nfaces'] = int(nfaces)
        if '- Turbulent Navier-Stokes' in line:
            nextline = next(lineiterator)
            itrb = nextline.split('(')[-1].split(')')[0]
            datum['turbulence'] = itrb
            datum['neq'] += teqtable[itrb]
        if 'Vibrational energy relaxation enabled' in line:
            datum['neq']+=1
        if 'Species list:' in line:
            species = line.strip().split()[3]
            ns = len(species.split(','))
            datum['ns'] = ns
            datum['neq'] += ns
        if 'Finite rate chemistry enabled' in line:
            datum['chemistry'] = True
        if 'PrgEnv-cray' in line:
            datum['machine'] = 'magnus'

        if line.startswith('iter '):
            header = line.split()
            next(lineiterator) # Skip the --- -----   ---- stuff
            break

    else:
        datum['iters'] = 0      # This one has no lines, set the start
        datum['start_iter'] = 2 # and stop indecies to return nothing
        datum['stop_iter'] =  1 # when slicing the main data list later
        return datum
    
    # We're now through the header, record anything that looks
    # like convergence info, by checking iteration number.
    try: numbers = next(lineiterator).split()
    except(StopIteration): 
        return datum

    iterno = int(numbers[0])
    start_time = float(numbers[header.index('time')])
    startiterno = iterno
    datum['start_iter'] = startiterno
    datum['start_time'] = start_time
    savediter = startiterno
    saved_final_time = start_time
    final_time = start_time

    for line in lineiterator:
        if line.startswith(str(iterno+1)+" "):
            final_time = float(line.split()[header.index('time')])
            iterno+=1
        if line.startswith('== Writing HDF5'):
            savediter = iterno
            saved_final_time = final_time
        if line.startswith('-- Time :'):
            line=line.split()
            datum['time']=float(line[-2])*unitconv[line[-1]]
            datum['hours']=datum['time']/3600.0
        if 'Walltime Used:' in line and 'time' not in datum.keys():
            line=line.split()
            h,m,s = list(map(float,line[-1].split(':')))
            datum['time'] = unitconv['secs']*s + unitconv['mins']*m + unitconv['hrs']*h
        if 'Resource Usage on' in line:
            datum['machine'] = 'raijin'

    # Have checked all lines, work out many iterations was done.
    datum['iters']=iterno-startiterno+1
    datum['stop_iter'] = savediter
    datum['final_iter'] = iterno
    datum['final_time'] = saved_final_time

    return datum    

def rm_slices(filename, rmkeys):
    """ Remove some slices from a file. I seem to have fixed this """
    #raise Exception("https://xkcd.com/1421/")
    f = h5py.File(filename)
    try:

        data = f['results/cells/list_1']
        initnumber = data.attrs['inum'][0]
        datanames = ['data_{}'.format(i) for i in range(1,initnumber+1)]
        datakeys = list(data.keys())

        if len(set(rmkeys))!=len(rmkeys):
            raise Exception("Duplicates in rmkeys")
        if len(datakeys)!=len(datanames):
            raise Exception("Length mismatch: keys={} names={}".format(len(datakeys), len(datanames)))
        if len(rmkeys)>initnumber:
            raise Exception("Too many rmkeys")
        for key in rmkeys:
            if key not in datakeys:
                raise Exception("Key missing!: {}".format(key))
        if len(rmkeys)==0:
            raise Exception("No datasets to remove!")

        print("Data sets to be removed: ")
        for key in rmkeys:
            time = data[key].attrs['time'][0]
            ncr  = data[key].attrs['ncr'][0]
            print(key, "time: {} ncr: {}".format(time,ncr))

        inputted = 0
        while True:
            val = input("Proceed? (y/n): ")
            if val=='y':
                # First remove the offending entries
                for name in rmkeys:
                    print("Removing", name)
                    idx = datanames.index(name)
                    del data[name]
                    del datanames[idx]

                # Now cleanup 
                for name  in datanames:
                    newname = name.replace('data','temp')
                    print("Renaming:    ", name, newname)
                    data[newname] = data[name]
                    del data[name]
                for i,name in enumerate(datanames):
                    newname = name.replace('data','temp')
                    rename = 'data_'+str(i+1)
                    print("Renumbering: ", newname, rename )
                    data[rename] = data[newname]
                    del data[newname]

                data.attrs['inum'] = initnumber-len(rmkeys)
                print("... Done")
                break
            elif val=='n':
                print('Aborting!')
                break
            else:
                print("Input not recognised")
            inputted+=1
            if inputted>3:
                print("Aborting!")
                break

    finally:
        f.close()
    return

def trim_slices(filename, ncr):
    """ Remove all slices from a file laid down after iteration ncr """
    f = h5py.File(filename)
    datalist = f['results/cells/list_1']

    try:
        initnumber = datalist.attrs['inum']
        rmnames = []
        print("Beginning with {} datasets".format(initnumber))

        for name, dataset in datalist.items():
            if dataset.attrs['ncr']>ncr:
                rmnames.append(name)

        print("Data sets to be removed: ")
        for i in rmnames: print(i)
        print("\nFor a total of {}".format(len(rmnames)))

        print("Removing...")
        for name in rmnames:
            del datalist[name]
        datalist.attrs['inum'] = initnumber-len(rmnames)
        print("... Done")

    finally:
        f.close()
    return

def remove_duplicate_slices(filename):
    """ Function to remove slices from h5 files that appear at duplicated times """
    times  = scrape_attrs(filename, 'results/cells/list_1', 'time')
    ntimes = len(times)
    # It's very important you don't sort the times here, because we are referring to them by index
    
    diffs = [[abs(i-j)/i for i in times] for j in times]
    similar = [[diff<1e-3 and i!=j for i,diff in enumerate(diffj)] for j,diffj in enumerate(diffs)]
    dupes = [[i for i in range(ntimes) if similarj[i]] for similarj in similar]

    ndupes = sum(1 for i in dupes if i!=[])
    print("Found duplicates: ", ndupes)
    if ndupes==0: return

    rmidxs = []
    for i,dupe in enumerate(dupes):
        if dupe!=[]:
            print("Slice {} with time {} matches:".format(i,times[i]))
            for d in dupe:
                print("   match: {} with time {}".format(d,times[d])) 
            keepme = all([i<d for d in dupe])
            if not keepme:
                rmidxs.append(i)
            print(" ")

    rm_slices(filename, rmidxs)
            
    return 

def readcsv(filename):
    """ Parse paraview csv files. """
    with open(filename) as fp:
        fiter = iter(fp)
        header = next(fiter).strip().split(',')
        header = [h.strip('"') for h in header]

        lines = [line.split(',') for line in fiter]
        cols = zip(*lines)
        cols = [np.array(list(map(float,col))) for col in cols]
        
        if len(cols)!=len(header): raise Exception

    return dict(zip(header,cols)) # Turn rows into columns

def humantime(time, precision=2):
    """ Turn time floats into human readable numbers """
    if abs(time)>1.0:
        humantime = str(round(time/1.0, precision)) + 's'
    if abs(time)>1e-3:
        humantime = str(round(time/1e-3, precision)) + 'ms'
    if abs(time)>1e-6:
        humantime = str(round(time/1e-6, precision)) + 'us'
    else:
        humantime = str(round(time/1e-9, precision)) + 'ns'
    return humantime 

STANDARD_HEADER='iter      residual         mass-bal         time            timestep    CFL'.split()

def make_line(header, line):
    numbers = []
    for val in line.split():
        try:
           numberval = int(val)
        except(ValueError):
           numberval = float(val)
        numbers.append(numberval)
    return dict(zip(header, numbers))

def cleanup(filename, header=None):
    if header==None: header = STANDARD_HEADER

    with open(filename) as fp:
        lines = [line.strip() for line in fp]

    cleanlines = []

    total = int(lines[-1].split()[0])
    print("Maximum steps: ", total)

    iters = np.zeros(total)
    data  = np.zeros((total, len(header)-1))
    overwritten = 0
    processed = 0

    for line in lines:
        linedict = make_line(header, line)
        iter = linedict['iter']
        if iters[iter-1]!=0: overwritten+=1

        iters[iter-1] = iter            # This is a bit silly.
        processed+=1

        for k,v in linedict.items():
            if k=='iter': continue
            data[iter-1,header.index(k)-1] = linedict[k]


    print("Done {} lines. {} overwritten...".format(processed, overwritten))
    return iters, data


if __name__=='__main__':
    #data = extract_data('us3dJul16-*.log')
    #print data.keys()
    #print len(data['iter'])
    #for k in data:
        #print data[k][:5]

    metadata=extract_metadata('us3dJul16-142*.log','us3dJul10-1111.log')
    for i in metadata:
        for k,v in i.items():
            print(k,v)


