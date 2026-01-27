import numpy as np 
from pathlib import Path

x=np.linspace(0,2048,64)
y=np.linspace(0,2048,64)
X,Y=np.meshgrid(x,y)
Z=np.zeros(X.shape)
target=open(Path("terrain.roughness.old").as_posix(),"w")
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        r=np.sqrt((X[i,j]-1024)**2+(Y[i,j]-1024)**2)
        if(r<=256):
            Z[i,j]=1e-4
        else:
            Z[i,j]=0.1
        target.write(f"{X[i,j]} {Y[i,j]} {Z[i,j]}\n")
target.close()

data = np.loadtxt(Path("terrain.roughness.old").as_posix())
x = np.unique(data[:, 0])
y = np.unique(data[:, 1])
z = data[:, 2]
assert len(z) == (len(x) * len(y))
with open(Path("terrain.roughness").as_posix(), "w") as f:
    f.write(f"{len(x)}\n")
    f.write(f"{len(y)}\n")
    x.tofile(f, sep="\n")
    f.write(f"\n")
    y.tofile(f, sep="\n")
    f.write(f"\n")
    z.tofile(f, sep="\n")

import matplotlib.pyplot as plt
plt.contourf(X, Y, Z, levels=50, cmap='viridis')
plt.colorbar(label='Roughness')
plt.show()
