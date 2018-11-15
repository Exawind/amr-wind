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
> cd exec/periodic_vortices
> make -j
> ./incflo3d.gnu.MPI.ex inputs
```

# Contributing

We welcome contributions in the form of pull-requests from anyone.  

# Acknowledgement of external libraries

This code relies on [Blitz++](https://github.com/blitzpp/blitz/wiki) and
[Algoim](https://fastmath-scidac.llnl.gov/software/algoim.html). 
