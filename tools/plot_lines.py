import numpy as np
from matplotlib import pylab
from sys import argv
import os
import shutil
import glob
import pandas as pd


# inputs 
# an average utau 
utau =  0.53
height = 1
vel = 8 
T0 = 300
kappa = 0.4
nondim = 1

a=pd.read_csv("averages.csv");

utau2 = utau*utau
#if nondim == 1:
a.z /= height
a.uu /= utau2
a.vv /= utau2
a.ww /= utau2
a.uv /= utau2
a.uw /= utau2
a.vw /= utau2


wind_dir =  180.0 + np.arctan2(a.u,a.v)*180.0/np.pi


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

umag = np.sqrt(a.u*a.u + a.v*a.v)
dudz = np.gradient(umag,a.z[1]-a.z[0])
#dudz = np.gradient(a.u,a.z[1]-a.z[0])
pylab.clf()
pylab.plot(a.z*kappa/utau*dudz,a.z, label=r"$\phi$")
pylab.ylim(0, 600.0/height)
pylab.xlim(0, 3)
pylab.xlabel(r"$\frac{\kappa z}{u_\tau} \frac{\partial u}{ \partial z}$")
pylab.ylabel("Height (m)")
#pylab.legend()
pylab.grid()
pylab.savefig("phi.pdf")


pylab.clf()
pylab.plot(a.wuu,a.z, label="<w'u'u'>")
pylab.plot(a.wuv,a.z, label="<w'u'v'>")
pylab.plot(a.wuw,a.z, label="<w'u'w'>")
pylab.plot(a.wvv,a.z, label="<w'v'v'>")
pylab.plot(a.wvw,a.z, label="<w'v'w'>")
pylab.plot(a.www,a.z, label="<w'w'w'>")
pylab.xlabel(r"$<u_iu_iu_i>$")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("www.pdf")


pylab.clf()
pylab.plot(a.www/np.float_power(a.ww,3.0/2.0),a.z, label="<w'w'w'>")
pylab.xlabel(r"$<u_iu_iu_i>/<u_iu_i>^{3/2}$")
pylab.ylabel("Height (m)")
#pylab.legend()
pylab.grid()
pylab.savefig("www_o_ww.pdf")


pylab.clf()
pylab.plot(a.Tu,a.z, label="<T'u'>")
pylab.plot(a.Tv,a.z, label="<T'v'>")
pylab.plot(a.Tw,a.z, label="<T'w'>")
pylab.xlabel(r"$<Tu_i>$")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("Tu.pdf")
