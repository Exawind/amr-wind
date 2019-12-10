import numpy as np
from matplotlib import pylab
from sys import argv
import os
import shutil
import glob
import pandas as pd


# inputs 
# an average utau 
utau = 0.5051855757
height = 1000
vel = 9
T0 = 300
nondim = 1

a=pd.read_csv("averages.csv");

utau2 = utau*utau
if nondim == 1:
    a.z /= height
    a.u /= vel
    a.v /= vel
    a.w /= vel
    a.uu /= utau2
    a.vv /= utau2
    a.ww /= utau2
    a.uv /= utau2
    a.uw /= utau2
    a.vw /= utau2
    a.theta /= T0

wind_dir =  360.0 + np.arctan(a.v/a.u)*180.0/np.pi


pylab.plot(np.sqrt(np.multiply(a.u,a.u) + np.multiply(a.v,a.v)),a.z, label='|U|')
pylab.plot(a.u,a.z, label=r'$U_x$')
pylab.plot(a.v,a.z, label=r'$U_y$')
pylab.xlabel("Velocity (m/s)")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("velmean.pdf")


pylab.clf()
pylab.plot(a.w,a.z, label=r'$U_z$')
pylab.xlabel("Velocity (m/s)")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("wmean.pdf")



pylab.clf()
pylab.plot(wind_dir,a.z, label='wind direction')
pylab.xlabel("Wind direction (deg)")
pylab.ylabel("Height (m)")
#pylab.legend()
pylab.grid()
pylab.savefig("wind_dir.pdf")




pylab.clf()
pylab.plot(a.uu,a.z, label="<u'u'>")
pylab.plot(a.vv,a.z, label="<v'v'>")
pylab.plot(a.ww,a.z, label="<w'w'>")
pylab.xlabel(r"$<u_iu_i>/u_\tau^2$")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("uu.pdf")


pylab.clf()
pylab.plot(a.uv,a.z, label="<u'v'>")
pylab.plot(a.uw,a.z, label="<u'w'>")
pylab.plot(a.vw,a.z, label="<v'w'>")
pylab.xlabel(r"$<u_iu_i>/u_\tau^2$")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("uv.pdf")


pylab.clf()
pylab.plot(a.theta,a.z, label="<T>")
pylab.xlabel("Temperature [K]")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("T.pdf")
