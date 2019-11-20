import yt
import numpy as np
from matplotlib import pylab

filename = "run128/line_plot"

# load one file to get the number of cells in 1D
ds = yt.load(filename+str(1).zfill(5))

# ds.print_stats()
# ds.field_list
# ds.derived_field_list
# ds.all_data()

g = ds.index.grids[0]
z = np.array(g["z"][0,0,:])
n = len(z)

u = np.zeros(n)
v = np.zeros(n)
w = np.zeros(n)
uu = np.zeros(n)
vv = np.zeros(n)
ww = np.zeros(n)
uv = np.zeros(n)
uw = np.zeros(n)
vw = np.zeros(n)
T = np.zeros(n)

start = 60000 
end = 82000
skip = 20
N = 0.0

for i in range(start,end+1,skip):
    ds = yt.load(filename+str(i).zfill(5));
    g = ds.index.grids[0]
    u += g["<u>"][0,0,:]
    v += g["<v>"][0,0,:]
    w += g["<w>"][0,0,:]
    uu += g["<u'u'>"][0,0,:]
    vv += g["<v'v'>"][0,0,:]
    ww += g["<w'w'>"][0,0,:]
    uv += g["<u'v'>"][0,0,:]
    uw += g["<u'w'>"][0,0,:]
    vw += g["<v'w'>"][0,0,:]
    T += g["<T>"][0,0,:]
    N += 1.0

# input an average utau 
utau = .25
utau2 = utau*utau

u /= N
v /= N
w /= N
uu /= N*utau2
vv /= N*utau2
ww /= N*utau2
uv /= N*utau2
uw /= N*utau2
vw /= N*utau2
T /= N

wind_dir =  360.0 + np.arctan(v/u)*180.0/np.pi


pylab.plot(np.sqrt(np.multiply(u,u) + np.multiply(v,v)),z, label='|U|')
pylab.plot(u,z, label=r'$U_x$')
pylab.plot(v,z, label=r'$U_y$')
pylab.xlabel("Velocity (m/s)")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("velmean.pdf")


pylab.clf()
pylab.plot(w,z, label=r'$U_z$')
pylab.xlabel("Velocity (m/s)")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("wmean.pdf")



pylab.clf()
pylab.plot(wind_dir,z, label='wind direction')
pylab.xlabel("Wind direction (deg)")
pylab.ylabel("Height (m)")
#pylab.legend()
pylab.grid()
pylab.savefig("wind_dir.pdf")




pylab.clf()
pylab.plot(uu,z, label="<u'u'>")
pylab.plot(vv,z, label="<v'v'>")
pylab.plot(ww,z, label="<w'w'>")
pylab.xlabel(r"$<u_iu_i>/u_\tau^2$")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("uu.pdf")


pylab.clf()
pylab.plot(uv,z, label="<u'v'>")
pylab.plot(uw,z, label="<u'w'>")
pylab.plot(vw,z, label="<v'w'>")
pylab.xlabel(r"$<u_iu_i>/u_\tau^2$")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("uv.pdf")


pylab.clf()
pylab.plot(T,z, label="<T>")
pylab.xlabel("Temperature [K]")
pylab.ylabel("Height (m)")
pylab.legend()
pylab.grid()
pylab.savefig("T.pdf")
