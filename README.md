# AMR-Wind 

[Website](https://www.exawind.org/) | [Documentation](https://amr-wind.readthedocs.io) | [Nightly test dashboard](http://my.cdash.org/index.php?project=AMR-Wind) 

[![Powered by AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io/amrex/) [![Build Status](https://github.com/Exawind/amr-wind/workflows/AMR-Wind-CI/badge.svg)](https://github.com/Exawind/amr-wind/actions) [![Docs Status](https://readthedocs.org/projects/pip/badge/?version=latest)](https://amr-wind.readthedocs.io)


AMR-Wind is a massively parallel, block-structured adaptive-mesh, incompressible
flow sover for wind turbine and wind farm simulations. The codebase is a
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

Documentation is available online at http://amr-wind.readthedocs.io/en/latest/.

## Compilation and usage

AMR-Wind is built upon the [AMReX library](https://amrex-codes.github.io/amrex).
A snapshot of the AMReX library is distributed along with the AMR-Wind source
code as a `git-submodule`. In addition to the AMReX library, you will require a
modern C++ compiler that supports the C++14 standard. Users wishing to execute
the code on high-performance computing (HPC) systems will also need MPI
libraries installed on their system. The code can also be compiled using NVIDIA
CUDA to target NVIDIA GPUs.

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

