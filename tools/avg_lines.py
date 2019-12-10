import yt
import numpy as np
from matplotlib import pylab
from sys import argv
import os
import shutil
import glob
import pandas as pd
from mpi4py import MPI
import copy


# inputs 
# an average utau 
start = 80000
end = 80010
skip = 5

filename = "run2304/line_plot"

fns = []
for i in range(start,end+1,skip):
    fn = filename+str(i).zfill(5);
    fns.append(fn)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

lfns = fns[rank::nprocs]
print(lfns)

# load one file to get the number of cells in 1D
ds = yt.load(fns[0])

# ds.print_stats()
# ds.field_list
# ds.derived_field_list
# ds.all_data()

g = ds.index.grids[0]
z = np.array(g["z"].flatten())
n = len(z)
a=pd.DataFrame(np.zeros((n,11)),columns=["z","u","v","w","uu","vv","ww","uv","uw","vw","theta"])
a.z += z[:]

print(a.z)
N = np.zeros(1)

for i, fn in enumerate(lfns):
    ds = yt.load(fn);
    g = ds.index.grids[0]
    a.u += g["<u>"].flatten() #[0,0,:]
    a.v += g["<v>"].flatten()
    a.w += g["<w>"].flatten()
    a.uu += g["<u'u'>"].flatten()
    a.vv += g["<v'v'>"].flatten()
    a.ww += g["<w'w'>"].flatten()
    a.uv += g["<u'v'>"].flatten()
    a.uw += g["<u'w'>"].flatten()
    a.vw += g["<v'w'>"].flatten()
    a.theta += g["<T>"].flatten()
    N[0] += 1.0

b = copy.copy(a.values)
# not sure how to do this with dataframe yet so for now use ndarray
comm.Reduce(a.values,b, MPI.SUM, 0)
a_avg = pd.DataFrame(b,columns=a.columns)

Ntotal = np.zeros(1)
comm.Reduce(N,Ntotal, MPI.SUM, 0)


if rank == 0:
    a_avg /= Ntotal[0]
    a_avg.z[:]= a.z[:]
    print(a_avg)
    a_avg.to_csv("averages.csv", index=False)
