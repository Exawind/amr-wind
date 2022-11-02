# AMR-Wind 

[Website](https://www.exawind.org/) | [User manual](https://exawind.github.io/amr-wind) | [API docs](https://exawind.github.io/amr-wind/api_docs) | [Nightly test dashboard](http://my.cdash.org/index.php?project=AMR-Wind) 

[![Powered by AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io/amrex/) [![Build Status](https://github.com/Exawind/amr-wind/workflows/AMR-Wind-CI/badge.svg)](https://github.com/Exawind/amr-wind/actions)


AMR-Wind is a massively parallel, block-structured adaptive-mesh, incompressible
flow solver for wind turbine and wind farm simulations. The codebase is a
wind-focused fork of [incflo](https://github.com/AMReX-Codes/incflo). The solver
is built on top of the [AMReX library](https://amrex-codes.github.io/amrex).
AMReX library provides the mesh data structures, mesh adaptivity, as well as the
linear solvers used for solving the governing equations. AMR-Wind is actively
developed and maintained by a dedicated multi-institutional team from [Lawrence
Berkeley National Laboratory](https://www.lbl.gov/), [National Renewable Energy
Laboratory](https://nrel.gov), and [Sandia National
Laboratories](https://sandia.gov).

The primary applications for AMR-Wind are: performing large-eddy simulations
(LES) of atmospheric boundary layer (ABL) flows, simulating wind farm
turbine-wake interactions using actuator disk or actuator line models for
turbines, and as a background solver when coupled with a near-body solver (e.g.,
[Nalu-Wind](https://github.com/exawind/nalu-wind)) with overset methodology to
perform blade-resolved simulations of multiple wind turbines within a wind farm.
For offshore applications, the ability to model the air-sea interaction effects
and its impact on the ABL characteristics is another focus for the code
development effort. As with other codes in the
[Exawind](https://github.com/exawind) ecosystem, AMR-wind shares the following
objectives:

- an open, well-documented implementation of the state-of-the-art computational
  models for modeling wind farm flow physics at various fidelities that are
  backed by a comprehensive verification and validation (V&V) process;

- be capable of performing the highest-fidelity simulations of flowfields within
  wind farms; and 

- be able to leverage the high-performance leadership class computating
  facilities available at DOE national laboratories.

## Documentation

Documentation is organized into a [user manual](https://exawind.github.io/amr-wind)
and a developer-focused [API
documentation](https://exawind.github.io/amr-wind). You can either
browse the docs online by following the links, or you can generate them locally
after downloading the code. Please follow the instructions in user manual to
build documentation locally.

## Compilation and usage

AMR-Wind is built upon the [AMReX library](https://amrex-codes.github.io/amrex).
A snapshot of the AMReX library is distributed along with the AMR-Wind source
code as a `git-submodule`. In addition to the AMReX library, you will require a
modern C++ compiler that supports the C++17 standard. Users wishing to execute
the code on high-performance computing (HPC) systems will also need MPI
libraries installed on their system. The code can also be compiled using MPI+X, 
where X can be OpenMP for CPU shared memory parallelism,
CUDA to target NVIDIA GPUs, HIP for AMD GPUs, or DPC++ for Intel GPUs.

### Contributing, reporting bugs, and requesting help

To report issues or bugs please [create a new
issue](https://github.com/Exawind/amr-wind/issues/new) on GitHub.

We welcome contributions from the community in form of bug fixes, feature
enhancements, documentation updates, etc. All contributions are processed
through pull-requests on GitHub.

## License

AMR-Wind is licensed under BSD 3-clause license. Please see the
[LICENSE](https://github.com/Exawind/amr-wind/blob/development/LICENSE) included in
the source code repository for more details.

