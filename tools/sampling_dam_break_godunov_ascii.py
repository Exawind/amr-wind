#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

loc_dir = "."
pp_dir = loc_dir + "/post_processing"

nt = 5 + 1
out_int = 30

plt.figure()
base_color = 0.8

for n in range(nt):
    x_oo=np.genfromtxt(pp_dir+"/sampling"+str(n*out_int).zfill(5)+".txt",delimiter=' ',skip_header=5)
    # Reorder array
    ind_0 = np.argsort(x_oo[:,0]) 
    x = np.zeros(np.shape(x_oo))
    flag = np.zeros(np.shape(x_oo)[0])
    for i in range(0, len(ind_0)): 
        x[i,:]= x_oo[ind_0[i],:]
        # When interface is not present, z location is set to zlo
        flag[i] = 1.0 if x[i,2] > 1e-8 else 0.
  
    # Shorten arrays (exclude points where interface not detected)
    nvalid = int(flag.sum())
    # Columns are position (3), id, cpu, ints of struct (uid, sid, nid), field components
    xshort = x[0:nvalid,0]
    zshort = x[0:nvalid,2]
    vshort = np.sqrt(x[0:nvalid,8]**2+x[0:nvalid,9]**2+x[0:nvalid,10]**2)

    color = base_color - (base_color) * (n/nt)
    cstr = str(color)

    plt.plot(xshort,zshort,color=cstr)
    plt.scatter(xshort,zshort,c=vshort,cmap="jet",vmin=0.,vmax=2.)

plt.ylabel('$z$',fontsize=16)
plt.xlabel('$x$',fontsize=16)
plt.colorbar()
plt.savefig('plot_sampling_ascii.pdf',format='pdf',dpi=300,bbox_inches="tight")