.. _capabilities:

Capabilities and Roadmap
========================

This section documents a non-exhaustive list of current AMR-Wind
capabilities and roadmap for future capabilities.

.. tip::

   If your project relies on a capability that is not yet present in
   AMR-Wind, please create an issue on the code project page.

Capabilities
------------

* Numerical methods

   * Advection: second order, piecewise parabolic, piecewise linear, WENO, Bell-Dawson-Shubin.

   * Diffusion: second order, explicit, Crank-Nicolson, and implicit

   * Mesh refinement: static mesh refinement (specified regions,
     geometric base refinement), adaptive mesh refinement (e.g., field
     based, curvature, q-criterion, vorticity)

   * Mesh mapping for non-uniform cartesian grids

* Equations systems

   * Incompressible and low Mach formulations of Navier-Stokes

   * Temperature

   * Level set

   * Subgrid scale kinetic energy

   * Specific dissipation rate

   * Passive scalar

* Turbulence modeling

   * Large Eddy Simulation: constant Smagorinsky,  AMD, one equation :math:`k_{sgs}`, non-linear, Kosovic

   * Wall models: log-law, constant stress, Schumann,

   * Reynolds-Average Navier-Stokes: :math:`k`-:math:`\omega` SST (and IDDES variant)

* Wind energy physics

   * Atmospheric boundary layer (ABL): various stability states
     (stable, unstable, neutral), precursor simulations with inflow
     boundary planes for wind farm simulations, anelastic formulation,
     mesoscale forcing, geostrophic forcing, Coriolis forcing,
     Monin-Obukhov similarity theory.

   * Actuator turbine representations: Joukowsky disks, uniform disks, actuator line

   * Coupling with OpenFAST

   * Coupling with Nalu-Wind for blade resolved simulations

* Postprocessing:

   * Visualization with VisIt, Paraview, yt

   * Sampling of fields with planes, point probes, lines, volumes, lidar, and radar

   * Scalar outputs such as kinetic energy, enstrophy, total wave energy, and norms

   * Free surface for gas-liquid interfaces

   * Turbulence averaging quantities such as Reynolds stresses

   * Field plane averaging and second and third order moments

   * Derived fields and field operators such as vorticity, q-criterion, strain-rates, gradients, divergence, laplacian

   * in-situ post-processing with Ascent

* Boundary conditions

   * Periodic, outflow, inflow, walls, user-defined inflows

   * Wall models (e.g., wall functions, stress)

   * Inflow planes from precursor simulations

   * Mesoscale forcing

* Multiphase flows

   * todo

* Geometry

   * Immersed boundary

   * Coupling with Nalu-Wind for body-conforming meshes with overset methodology

* Transport models

   * Constant transport coefficients

   * Two phase transport

* High performance computing

   * Shared memory parallelism with OpenMP threading

   * Distributed memory parallelism with MPI

   * Supports all major compilers (e.g., GCC, Intel, LLVM)

   * Runs on all major GPU vendors (NVIDIA, AMD, Intel)

   * Supported build systems: cmake, spack

* Linear solvers

   * native AMReX solvers such as MLMG

   * hypre

Roadmap
-------

The roadmap is an evolving, living document and does not purport to
track every future capability. It is not a promise of future
capabilities. The main use case is to inform users of
potential upcoming new capabilities.
