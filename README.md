# Directory overview

| File       | Description                                         |
| -----------| --------------------------------------------------- |
| test       | Directory for building EB executable                |
| test_no_eb | Directory for building non-EB executable            |
| src        | C++/Fortran source files                            |


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
> cd test
> make -j4
> mpirun -np 4 incflo3d.gnu.MPI.ex inputs.channel_cylinder
```

# Contributing

We welcome contributions in the form of pull-requests from anyone.
