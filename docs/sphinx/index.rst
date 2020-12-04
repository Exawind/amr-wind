==========
 AMR-Wind
==========

`AMR-Wind <https://github.com/exawind/amr-wind>`_ is a massively parallel,
block-structured adaptive-mesh, incompressible flow sover for wind turbine and
wind farm simulations. The codebase is a wind-focused fork of `incflo
<https://github.com/AMReX-Codes/incflo>`_. The solver is built on top of the
`AMReX library <https://amrex-codes.github.io/amrex>`_. AMReX library provides
the mesh data structures, mesh adaptivity, as well as the linear solvers used
for solving the governing equations. AMR-Wind is actively developed and
maintained by a dedicated multi-institutional team from `Lawrence Berkeley
National Laboratory <https://www.lbl.gov/>`_, `National Renewable Energy
Laboratory <https://nrel.gov>`_, and `Sandia National Laboratories
<https://sandia.gov>`_.

The primary applications for AMR-Wind are: performing large-eddy simulations
(LES) of atmospheric boundary layer (ABL) flows, simulating wind farm
turbine-wake interactions using actuator disk or actuator line models for
turbines, and as a background solver when coupled with a near-body solver (e.g.,
`Nalu-Wind <https://github.com/exawind/nalu-wind>`_) with overset methodology to
perform blade-resolved simulations of multiple wind turbines within a wind farm.
For offshore applications, the ability to model the air-sea interaction effects
and its impact on the ABL characteristics is another focus for the code
development effort. As with other codes in the
`Exawind <https://github.com/exawind>`_ ecosystem, AMR-wind shares the following
objectives:

- an open, well-documented implementation of the state-of-the-art computational
  models for modeling wind farm flow physics at various fidelities that are
  backed by a comprehensive verification and validation (V&V) process;

- be capable of performing the highest-fidelity simulations of flowfields within
  wind farms; and

- be able to leverage the high-performance leadership class computating
  facilities available at DOE national laboratories.


.. toctree::
   :maxdepth: 2

   user/user
   theory/theory 
   developer/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

