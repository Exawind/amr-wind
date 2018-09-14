# Directory overview

| File      | Description                                         |
| ----------| --------------------------------------------------- |
| doc       | Documentation                                       |
| exec      | Run directory for executables                       |
| src       | C++/Fortran source files                            |
| tools     | CMake configuration files                           |


# Using incflo

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
