#!/usr/bin/env python3

import numpy as np
import netCDF4 as ncdf
import matplotlib.pyplot as plt

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

loc_dir = "."
pp_dir = loc_dir + "/post_processing"

nt = 5 + 1
out_int = 30

class SamplingFile(object):
    """Interface to Sampling NetCDF file"""

    def __init__(self, sampling_file = "sampling.nc", group_name = "s1"):
        """
        Args:
            stats_file (path): Absolute path to the NetCDF file
        """
        self.sampling_file = sampling_file
        self.sampling = ncdf.Dataset(self.sampling_file)
        self.fs = self.sampling["/"+group_name]
        self._points = self.fs.variables["points"][:,:,:]
        self._velocityx = self.fs.variables["velocityx"][:,:]
        self._velocityy = self.fs.variables["velocityy"][:,:]
        self._velocityz = self.fs.variables["velocityz"][:,:]

    @property
    def locations(self):
        return self._points

    @property
    def vmag(self):
        return np.sqrt(self._velocityx**2+self._velocityy**2+self._velocityz**2)

plt.figure()
base_color = 0.8

sampl = SamplingFile(sampling_file=pp_dir+"/sampling00000.nc", \
                    group_name = "fs")

for n in range(nt):
    x_oo=sampl.locations[n,:,:]
    Vmag_oo = sampl.vmag[n,:]

    # Reorder arrays
    ind_0 = np.argsort(x_oo[:,0]) 
    x = np.zeros(np.shape(x_oo)[0])
    z = np.zeros(np.shape(x_oo)[0])
    Vmag = np.zeros(np.shape(x_oo)[0])
    flag = np.zeros(np.shape(x_oo)[0])
    for i in range(0, len(ind_0)): 
        x[i]= x_oo[ind_0[i],0]
        z[i]= x_oo[ind_0[i],2]
        Vmag[i] = Vmag_oo[ind_0[i]]
        # When interface is not present, z location is set to 0
        flag[i] = 1.0 if z[i] > 1e-8 else 0.
    # Shorten arrays (exclude points where interface not detected)
    nvalid = int(flag.sum())
    xshort = x[0:nvalid]
    zshort = z[0:nvalid]
    vshort = Vmag[0:nvalid]

    color = base_color - (base_color) * (n/nt)
    cstr = str(color)

    plt.plot(xshort,zshort,color=cstr)
    plt.scatter(xshort,zshort,c=vshort,cmap="jet",vmin=0.,vmax=2.)

plt.ylabel(r'$z$',fontsize=16)
plt.xlabel(r'$x$',fontsize=16)
plt.colorbar()
plt.savefig('plot_sampling_netcdf.pdf',format='pdf',dpi=300,bbox_inches="tight")