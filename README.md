# Directory overview

| File      | Description                                         |
| ----------| --------------------------------------------------- |
| exec      | Run directory for executables                       |
| src       | C++/Fortran source files                            |


# Using incflo

## Install [Blitz++](https://github.com/blitzpp/blitz/wiki)

APT (Debian):
```shell
> sudo apt install libblitz-doc libblitz0-dev libblitz0ldbl
```

[Homebrew](https://brew.sh/) (macOS):
```shell
> brew install blitz
```

## Build AMReX Library

Clone AMReX from the official Git repository and checkout the _development_ branch.
```shell
> git clone https://github.com/AMReX-Codes/amrex.git
> cd amrex
> git checkout development
```

## Build and run an example incflo problem
Clone and build incflo
```shell
> git clone http://github.com/AMReX-Codes/incflo.git
> cd exec
> make -j4
> mpirun -np 4 incflo3d.gnu.MPI.ex inputs.channel_cylinder
```

# Contributing

We welcome contributions in the form of pull-requests from anyone.  

# Acknowledgement of external libraries

Apart from 
[AMReX](https://amrex-codes.github.io/), this code relies on
[Blitz++](https://github.com/blitzpp/blitz/wiki) and
[Algoim](https://fastmath-scidac.llnl.gov/software/algoim.html). 
