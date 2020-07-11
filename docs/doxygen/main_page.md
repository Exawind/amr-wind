# AMR-Wind API documentation {#mainpage}

[AMR-Wind](https://github.com/exawind/amr-wind) is a massively parallel,
block-structured adaptive-mesh, incompressible flow sover for wind turbine and
wind farm simulations. The codebase is a wind-focused fork of
[incflo](https://github.com/AMReX-Codes/incflo). The solver is built on top of
the [AMReX library](https://amrex-codes.github.io/amrex). AMReX library provides
the mesh data structures, mesh adaptivity, as well as the linear solvers used
for solving the governing equations. AMR-Wind is actively developed and
maintained by a dedicated multi-institutional team from [Lawrence Berkeley
National Laboratory](https://www.lbl.gov/), [National Renewable Energy
Laboratory](https://nrel.gov), and [Sandia National
Laboratories](https://sandia.gov).

This document is intended for developers who want to understand the C++ code
structure and modify the codebase and, therefore, assumes that the reader is
familiar with the installation, compilation, and execution steps. If you are new to
AMR-Wind and haven't installed/used AMR-Wind previously, we recommend starting
with the [user manual](https://amr-wind.readthedocs.io) that provides a detailed
overview of the installation process as well as general usage.

## How to use this API guide?

This section provides a brief overview of the organization of this API guide so
as to enable readers to quickly find the sections that they are interested in.
The source code documentation is organized in [sections](modules.html) that
divides the codebase into logical groups. We recommend that you start
[here](modules.html) and navigate to the sections that you are interested in.

AMR-Wind is built on top of the [AMReX
library](https://amrex-codes.github.io/amrex/). The \ref core "core data structures" 
provide higher-level abstractions on top of AMReX data structures.
We recommend reading the [AMReX basics
chapter](https://amrex-codes.github.io/amrex/docs_html/Basics.html) to
familiarize yourself with the core AMReX terminology and concepts. Once you have
read that chapter, read the \ref core "AMR-Wind core"
documentation and familiarize yourself with the concept of \ref
amr_wind::Field "Field" and \ref amr_wind::FieldRepo "FieldRepo" (see \ref
fields) in AMR-Wind as these are used quite heavily everywhere in the code. Two
other global data structures that are used frequently are \ref amr_wind::CFDSim
"CFDSim" and \ref amr_wind::SimTime "SimTime". `CFDSim` represents the
simulation environment and holds references to the mesh, the field repository,
time instance, the registered \ref physics "physics" instances, the \ref eqsys
"equation systems", and \ref utilities "post-processing and I/O" utilities.
`SimTime` holds all attributes related to time within the code and determines
when to advance the simulation, exit, or write outputs.

### Source code organization

Upon successful download/clone, the base repository (`amr-wind`) has source code
is organized in subdirectories described below:

- `amr-wind` -- C++ source files. All code is located within this directory
- `unit_tests` -- Unit-tests for individual modules/classes
- `cmake` -- Functions, utilities used during CMake configuration phase
- `docs` -- User manual (Sphinx-based) and Doxygen files
- `submods` -- Third-party libraries and dependencies
- `test` -- Regression tests and associated infrastructure
- `tools` -- Miscellanous post-processing scripts and other utilities

When developing new features, we strongly recommend creating a unit-test and
develop features incrementally and testing as you add capabilities. Unit-tests
are also a good way to explore the usage of individual components of the code.

### Building API documentation locally

The API documentation is automatically generated from specially-formatted
comments in the source code using [Doxygen](https://www.doxygen.nl/index.html).
Please consult the [Doxygen
manual](https://www.doxygen.nl/manual/docblocks.html) to learn about documenting
code. To generate this documentation on your local machine, or to rebuild docs
during code development process you'll need to install `doxygen` and `graphviz`
executables on your system. Once these executables have been successfully
installed on your system, you can generate this documentation by executing the
following commands:

~~~~~~~~~~~.sh
git clone https://github.com/exawind/amr-wind.git
cd amr-wind

# Doxygen command should be executed from the top-level directory
doxygen docs/doxygen/Doxyfile
# The resulting documentation is in `build/html` directory

# Open main page on your browser
open build/html/index.html 
~~~~~~~~~~~
 
## Contributing

AMR-Wind is an open-source code and we welcome contributions from the community.
Please consult the [developer
documentation](https://amr-wind.readthedocs.io/en/latest/dev/index.html) section
of the user manual to learn about the process of submitting code enhancements,
bug-fixes, documentation updates, etc.

## License

AMR-Wind is licensed under BSD 3-clause license. Please see the
[LICENSE](https://github.com/Exawind/amr-wind/blob/development/LICENSE) included in
the source code repository for more details.

