#!/usr/bin/env python
from pylab import *

# ----
# Time averages the line_plot data given by dat, between times t1 and
# t2.
#
def timeaverage(dat, t1, t2, Ncells, Ncomp):
    t = dat[:,1]
    datfiltered=dat[(t>=t1)&(t<=t2),:]
    tstart = datfiltered[0,1]
    tend   = datfiltered[-1,1]
    istart = int(datfiltered[0,0])
    iend   = int(datfiltered[-1,0])
    Nlines = len(datfiltered)
    Nblocks = int(Nlines/Ncells)
    if (Nlines%Ncells) != 0:
        print("Something wrong with number of blocks")
    avgdat=zeros((Ncells, Ncomp))
    for i in range(Nblocks-1):
        i1 = i*Ncells
        i2 = (i+1)*Ncells
        i3 = (i+2)*Ncells
        dt = datfiltered[i2, 1] - datfiltered[i1, 1]
        avgdat=avgdat+0.5*dt*(datfiltered[i1:i2,:] + datfiltered[i2:i3,:])
    return avgdat/(tend-tstart)

# ----
# Loads the line_plot.txt file given by filename, for all times
# between t1 and t2.
#
# If the line_plot.txt file is missing a header, provde the input
# headerinfo=[ncell ncomp]
def loadfile(filename, t1, t2, headerinfo=[]):
    if len(headerinfo)>0: hasheader=False
    else: hasheader=True
    # Load the header info
    if hasheader:
        with open(filename) as fp:
            header = fp.readline()
            line   = fp.readline().strip().split(',')
            ncell  = int(line[0])
            ncomp  = int(line[1])
    else:
        ncell=headerinfo[0]
        ncomp=headerinfo[1]
    dat=[]
    with open(filename) as fp:
        if hasheader:  # Read past the headers
            junk      = fp.readline()
            junk      = fp.readline()
            headerline= fp.readline()
        colheaders=headerline.replace('#','',1).replace('\n','').split(',')
        # Read the result 
        line = fp.readline()
        iline= 1
        while line:
            splitline=line.strip().split(',')
            time=float(splitline[1])
            if (t1 <= time) and (time <= t2):
                dat.append([float(x) for x in splitline])
                #print("Adding time: %f"%time)
            if (time>t2):
                line=False
            else:
                line=fp.readline()
                iline=iline+1
    return array(dat), ncell, ncomp, colheaders

# ----
# Load the line_plot.txt file in filename, and average it from time t1
# to t2.  If line_plot.txt has no headers, supply headerinfo=[ncell
# ncomp]
#
def avglineplot(filename, t1,t2, headerinfo=[]):
    dat, ncell, ncomp, colheaders=loadfile(filename, t1, t2, headerinfo)
    tavgdat=timeaverage(dat, t1,t2, ncell, ncomp)
    return tavgdat, colheaders

# ----
# Returns the variable(s) matching the variables in varnames 
# Behavior depends on varnames:
# - If varnames is a single string (like varnames='z'), then returns
#   only that variable on output.
# - If varnames is a list of strings (varnames=['z', 'u_avg',
#   'v_avg']), then returns a dictionary with those variables.
#
def extractvars(dat, colheaders, varnames):
    if isinstance(varnames, list):   
        vardict=dict()
        for var in varnames:
            vardict[var] = dat[:,colheaders.index(var)]
        return vardict
    elif isinstance(varnames, str):
        return dat[:, colheaders.index(varnames)]
    else:
        print("varnames = "+repr(varnames)+" is not valid input")
    
# Split the data into z, u, v, w, and uu, vv, www columns
def splitdat(dat):
    z    = dat[:,2]
    u    = dat[:,3]
    v    = dat[:,4] 
    w    = dat[:,5]
    uu   = dat[:,8]
    vv   = dat[:,11] 
    ww   = dat[:,13]
    return z, u, v, w, uu, vv, ww

