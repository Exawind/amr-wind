import numpy as np 
from pathlib import Path


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