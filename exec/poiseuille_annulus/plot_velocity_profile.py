from matplotlib import pyplot as plt
import numpy as np
from platform import system

deltap = 6.7745e-3
mu = 1.0e-3
L = 0.5
Ro = 1
Ri = 0.25
alpha = Ri / Ro
C = deltap * Ro**2 / (4 * mu * L)
u = lambda r: C * (1 - (r / Ro)**2 + (1 - alpha**2) * np.log(Ro / r) / np.log(alpha))

# Set up figure environment
fig = plt.figure(num=1, figsize=(6,3))
fig.clf()

# x-velocity
datatype = np.dtype([('x', float), ('u', float)])
ax = fig.add_subplot(1,2,1)  
filename = 'visit.curve'
data = np.loadtxt(filename, datatype, skiprows=1)
dx = (data['x'][1] - data['x'][0])
x = data['x'] + dx / 2 - Ro
mask = (x > Ri)
ax.plot(x[mask], data['u'][mask], '.', markerfacecolor='none', label='Computed')
r = np.linspace(Ri, Ro, 1000)
ax.plot(r, u(r), label='Analytical')
ax.legend(loc='best')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')

# error
ax = fig.add_subplot(1,2,2)  
err = np.abs(u(x[mask]) - data['u'][mask])
relerr = err / u(x[mask])
ax.plot(x[mask], err, 'x', markerfacecolor='none', label='Absolute')
#ax.plot(x[mask], relerr, '+', markerfacecolor='none', label='Relative')
#ax.legend(loc='best')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'Error')

## Finalize
fig.tight_layout()
#plt.show()
plt.savefig("annulus.jpeg")
