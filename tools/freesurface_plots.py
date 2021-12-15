# -*- coding: utf-8 -*-
# pylint: disable=too-many-locals

"""
ABL Statistics plotting utility
"""

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import netCDF4 as ncdf
import pandas as pd
import sys


class FreeSurfaceStatsFile(object):
    """Interface to ABL Statistics NetCDF file"""

    def __init__(self, stats_file="FreeSurface.nc"):
        """
        Args:
            stats_file (path): Absolute path to the NetCDF file
        """
        self.stats_file = stats_file
        self.fs_stats = ncdf.Dataset(stats_file)
        self._time = self.fs_stats["time"][:]
        self._coords = self.fs_stats["coordinates2D"][:]
        self._heights = self.fs_stats["heights"][:]
        self.start = self.fs_stats.getncattr("start")
        self.end = self.fs_stats.getncattr("end")
        self.ninst = len(self.heights[0,:,0])
        ijk_dims = self.fs_stats.getncattr("ijk_dims")
        self.ni = ijk_dims[0]
        self.nj = ijk_dims[1]

    @property
    def time(self):
        """The time array from the FreeSurface file"""
        return self._time

    @property
    def coords(self):
        """The coordinates array where data is available"""
        return self._coords

    @property
    def heights(self):
        """Heights where data is available"""
        return self._heights

def plot_interface(stats, pdffile, num_step):
    """Generate interface plot"""
    ht_full = stats.heights
    coords = stats.coords
    time = stats.time
    start = stats.start
    end = stats.end
    nx = stats.ni
    ny = stats.nj
    ninst = stats.ninst
    print("Instant displayed = %.3f s"%(time[num_step]))
    ht_now = ht_full[num_step, :, :]

    # Plot 2D color plot for each instance that is available
    for n in range(ninst):
        plt.figure()
        arr = np.reshape(ht_now[n,:],(nx, ny))
        plt.imshow(arr, extent=[start[0],end[0],start[1],end[1]],
            origin="lower",interpolation='bilinear')
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.colorbar()
        pdffile.savefig(dpi=300, bbox_inches="tight")
        plt.close()

if __name__ == "__main__":
    clarg = sys.argv
    if (len(clarg) != 2):
        print("Program cancelled. Please input desired timestep for plotting.")
        quit()
    num_step = int(sys.argv[1])
    print("Timestep: %d"%(num_step))
    statsfile = FreeSurfaceStatsFile()
    with PdfPages("interface_test.pdf") as pdf:
        plot_interface(statsfile, pdf, num_step)
