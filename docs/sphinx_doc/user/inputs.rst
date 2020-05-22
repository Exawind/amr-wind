.. _amrwind-inputs:

AMR-Wind inputs file
=====================

To run :program:`amr_wind`, the user must provide a text file containing inputs
describing the problem and any additional command-line arguments that override
the parameters in the input file for that particular invocation of the
executable.

.. code-block:: console

   # Parse input parameters from `inputs.abl` but change max_step to 20
   $ ./amr_wind inputs.abl time.max_step=20


The input file is a simple text file containing ``key = value`` entries for the
input parameters. The text file can include comments, any text beginning with
``#`` till the end of line (EOF) is interpreted as comments and ignored by the
parser. Input file processing is handled by `AMReX ParmParse library
<https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`_. This
section documents the various input file parameters and their default values (if
available). In :program:`amr_wind`, the input file is broken into *sections*
indicated by a namespace prefix. For example, all inputs related to the problem
domain are prefixed with ``geometry.`` and so on. A sample input file is shown below

.. literalinclude:: ./amr_wind_inputs.txt
   :linenos:
   :emphasize-lines: 6, 11, 15, 26, 33, 49, 55, 59, 75, 85, 89, 95, 107, 111


Input file reference
---------------------

The AMR-Wind input file is organized in the following sections

======================= ============================================================
Section                 Description
======================= ============================================================
``geometry``            Computational domain information
``amr``                 Mesh refinement controls
``time``                Simulation time controls
``io``                  Input/Output controls
``incflo``              CFD algorithm and physics controls
``transport``           Transport equation controls
``turbulence``          Turbulence model controls 
``ABL``                 Atmospheric boundary layer (ABL) controls
``Momentum sources``    Activate Momentum source terms and their parameters
``Static Refinement``   Cartesian static refinement options
``Boundary conditions`` Boundary condition types and gradients
``MLMG options``        Multi-Level Multi-Grid Linear solver options
======================= ============================================================

This section documents the parameters available within each section.

.. note::

   - Boolean flags (``true/false``) can also be indicated using integers in the text file
     and uses the convention ``0 = False`` and ``1 = True``.

   - Quotes around strings are optional.

   - If an input is repeated only the last one is used.

Section: ``geometry``
~~~~~~~~~~~~~~~~~~~~~

This section deals with inputs related to the problem domain.

.. input_param:: geometry.prob_lo

   **type:** List of 3 real numbers, mandatory

   The coordinates of *lower corner* of the computational domain bounding box.

.. input_param:: geometry.prob_hi

   **type:** List of 3 real numbers, mandatory

   The coordinates of the *upper corner* of the computational domain bounding box.

.. input_param:: geometry.is_periodic

   **type:** List of 3 integers, mandatory

   Flags indicating whether the flow is periodic in the ``x``, ``y``, or ``z``
   directions respectively.

Section: ``amr``
~~~~~~~~~~~~~~~~

This section mostly contains inputs the adaptive mesh refinement capabilities of
:program:`amr_wind`.

.. input_param:: amr.n_cell

   **type:** List of 3 integers, mandatory

   The number of cells in the coarset level of AMR hierarchy

.. input_param:: amr.max_level

   **type:** Integer, optional, default: 0

   The maximum AMR level in the refinement hierarchy. Default value is ``0``
   indicating a single mesh level with uniform resolution in the three
   directions.

.. input_param:: amr.max_grid_size

   **type:** Integer, optional, default: 32

   Maximum number of cells at level 0 in each grid.
   There are also options to specify this value in each direction,
   please refer to AMReX documentation.

.. input_param:: amr.blocking_factor

   **type:** Integer, optional, default: 8

   Each grid must be divisible by ``amr.blocking_factor``.
   There are also options to specify this value in each direction,
   please refer to AMReX documentation.

Section: ``time``
~~~~~~~~~~~~~~~~~

This section deals with parameters that control the simulation.

.. input_param:: time.stop_time

   **type:** Real number

   The time (in seconds) when the simulation should stop. This parameter is
   related to :input_param:`time.max_step` depending on whether the user
   provides a positive or negative value for this input. If the user provides a
   positive value, then the stop criteria is such that the simulation stops when
   it reaches the ``stop_time``. If the user specifies a negative value, then
   the simulation stops when it has completed the prescribed number of steps as
   determined by :input_param:`time.max_step`.

.. input_param:: time.max_step

   **type:** Integer

   The number of timesteps to run before termination. See also
   :input_param:`time.stop_time`.

.. input_param:: time.fixed_dt

   **type:** Real number

   Fixed timestep size (in seconds) used to advance the simulation. If this
   parameter is negative, then :input_param:`time.cfl` is used to determine the
   adaptive timestep during the simulation.
   
.. input_param:: time.initial_dt

   **type:** Real number

   Initial timestep size (in seconds) used to initialize the simulation. 
   Only activated if :input_param:`time.fixed_dt` is negative 
   (signalling CFL controlled time stepping) and if ``time.initial_dt`` is positive.
   This parameter can be useful for starting CFL-controlled simulations like 
   Rayleigh-Taylor flow that initialize with zero velocity.

.. input_param:: time.cfl

   **type:** Real number, optional, default = 0.5

   The maximum CFL allowed during the simulation. Adaptive timestepping algorithm is 
   used to ensure that the maximum CFL condition is not violated while using the largest
   allowable timestep to advance the simulation.

.. input_param:: time.init_shrink

   **type:** Real number, optional, default = 0.1

   ``init_shrink`` indicates a reduction factor applied to the timestep size
   during initialization phased used to perform the initial iterations to
   determine pressure. This factor must be a value greater than zero but less
   than or equal to 1.0.

.. input_param:: time.regrid_interval

   **type:** Integer

   If :input_param:`amr.max_level` is greater than zero, this parameter
   indicates the frequency (in timesteps) at which the mesh is adaptively
   refined based on various user-specified criteria. If this value is negative,
   the mesh is only refined once during initialization and remains constant for
   the rest of the simulation.

.. input_param:: time.plot_interval

   **type:** Integer

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which outputs (plot files) are written to disk.

.. input_param:: time.checkpoint_interval

   **type:** Integer

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which checkpoint (restart) files are written to disk.

.. input_param:: time.use_force_cfl

   **type:** Boolean, optional, default = true

   If this flag is true then the forces (including the pressure gradient) are included
   in the CFL calculation.
   
Section: ``io``
~~~~~~~~~~~~~~~~~

This section deals with parameters that control input/output to the simulation.

.. input_param:: io.line_plot_int

   **type:** Integer, optional, default = 1

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which line plot files are written to disk. 
   Currently line plotting is only turned on for ABL simulations.
   
.. input_param:: io.line_plot_type

   **type:** Integer, optional, default = 0

   This input specifies the type of line plot output, where 0 outputs averages and
   fluctuations in a text format, 
   1 outputs only averages in text format, 
   and 2 outputs averages and fluctuations in a binary format. 
   Currently line plotting is only turned on for ABL simulations

   
.. input_param:: io.KE_int

   **type:** Integer, optional, default = -1

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which Kinetic Energy is computed and printed to the log file.
   
.. input_param:: io.check_file

   **type:** String, optional, default = "chk"

   If :input_param:`time.checkpoint_interval` is greater than zero this is the name of the checkpoint 
   file appended with the current timestep
   
.. input_param:: io.plot_file

   **type:** String, optional, default = "plt"

   If :input_param:`time.plot_interval` is greater than zero this is the name of the plot
   file appended with the current timestep
   
.. input_param:: io.restart_file

   **type:** String, optional, default = ""

   If a string is present `amr-wind` will restart using the specified file in the string.
   
   
Section: ``incflo``
~~~~~~~~~~~~~~~~~~~

This section deals with parameters that mostly determine how amr-wind is run such 
as initial conditions and discretization options.

.. input_param:: incflo.physics

   **type:** List of strings

   Specify a string or a list of strings for each type of physics to initialize and simulate.
   Physics is additive and more than one type of physics may be used.
   Current implemented physics are FreeStream, ABL, RayleighTaylor, BoussinesqBubble, 
   and TaylorGreenVortex
   
.. input_param:: incflo.density

   **type:** Real

   Specify reference density. 
   For the most part if :input_param:`incflo.constant_density` = true then `incflo.density` sets a constant density everywhere.
   Refer to the field initializer for your chosen :input_param:`incflo.physics` for how `incflo.density` is used.
   
.. input_param:: incflo.velocity

   **type:** List of three Real numbers, sometimes mandatory depending on the :input_param:`incflo.physics`

   Specify reference velocity in the "x-", "y-", "z-direction". 
   Refer to the field initializer for your chosen :input_param:`incflo.physics` for how `incflo.velocity` is used.
   
.. input_param:: incflo.delp

   **type:** List of 3 real numbers

   Specify a constant pressure change in the x,y, and z directions. 
   Currently this is broken and will abort if a non-zero number is specified.
   A user defined body force will replace this in the future.
   
.. input_param:: incflo.verbose

   **type:** Integer, optional, default = 0

   Specifies amount of verbosity. A value of 0 is minimal verbosity output and 3 gives full verbosity output. 
   
.. input_param:: incflo.initial_iterations

   **type:** Integer, optional, default = 3

   Number of initial pressure iterations to perform when the simulation starts.
   
.. input_param:: incflo.do_initial_proj

   **type:** Boolean, optional, default = true

   This flag when true performs a nodal projection
   to ensure the initial velocity is divergence-free. 
   
.. input_param:: incflo.constant_density

   **type:** Boolean, optional, default = true

   This flag specifies if density is constant throughout the simulation. 
   If the flag is true then a constant density field is used, the density field is copied from old to new time steps. 
   If the flag is false then density is not constant and an extra advection equation is solved to evolve density.
   
.. input_param:: incflo.use_godunov

   **type:** Boolean, optional, default = false

   Specifies which advection scheme to use: either method of lines (false) or Godunov (true). 
   The method of lines is the default option but Godunov is more accurate, 
   can handle a larger CFL number, and more computational efficient.
   
.. input_param:: incflo.use_ppm

   **type:** Boolean, optional, default = true

   When estimating the two states in a Godunov scheme a piecewise parabolic method (PPM) is used when this flag is true
   or when the flag is false a less accurate piecewise linear method (PLM) is used instead.
   Note: only used when :input_param:`incflo.use_godunov` = true.
   
.. input_param:: incflo.godunov_use_forces_in_trans

   **type:** Boolean, optional, default = false

   Specifies if body forces are included in the transverse velocity prediction.
   Note: only used when :input_param:`incflo.use_godunov` = true.
   
.. input_param:: incflo.diffusion_type

   **type:** Integer, optional, default = 2

   Determines how the diffusion term is handled when updating the momentum equations. 
   A value of 0 is explicit diffusion and all diffusion terms are moved to the right hand side 
   (warning this carries with it a more stringent CFL restriction), 
   a value of 1 is Crank-Nicolson and diffusion terms are on both the left and right hand sides,
   and a value of 2 (default) is a fully implicit diffusion where the entire diffusion term is handled on the left hand side.
   
.. input_param:: incflo.rhoerr

   **type:** Real number or a list of Real numbers

   When :input_param:`amr.max_level` > 0 this will trigger mesh adaption for density that is greater than `incflo.rhoerr`.
   This maybe specified as a single number for all levels or a value per AMR level.
   
.. input_param:: incflo.gradrhoerr

   **type:** Real number or a list of Real numbers

   When :input_param:`amr.max_level` > 0 this will trigger mesh adaption if the difference
   between density at a cell center and its neighbors is greater than `incflo.gradrhoerr`. 
   This maybe specified as a single number for all levels or a value per AMR level.
   
Section: ``transport``
~~~~~~~~~~~~~~~~~~~~~~

This section is for setting thermal and momentum diffusivities.

.. input_param:: transport.model

   **type:** String, optional, default = ConstTransport

   Currently this is the only transport model implemented.
   
.. input_param:: transport.viscosity

   **type:** Real, optional, default = 1.0e-5

   Sets the dynamic viscosity
   
.. input_param:: transport.laminar_prandtl 

   **type:** Real, optional, default = 1.0

   Sets the laminar Prandtl number.
   
.. input_param:: transport.turbulent_prandtl 

   **type:** Real, optional, default = 1.0

   Sets the turbulent Prandtl number.
   
Section: ``turbulence``
~~~~~~~~~~~~~~~~~~~~~~~

This section is for setting turbulence model parameters

.. input_param:: turbulence.model

   **type:** String, optional, default = Laminar

   Specifies which turbulence model to use, by default "Laminar" is chosen 
   (effectively no turbulence model). 
   Currently the only turbulence model is "Smagorinsky"

   
.. input_param:: Smagorinsky_coeffs.Cs

   **type:** Real, optional, default = 0.135

   Specifies the coefficient used in the `Smagorinsky` turbulence model. 
   
Section: ``ABL``
~~~~~~~~~~~~~~~~~~~~~~~

This section is for setting atmospheric boundary layer parameters.

.. input_param:: ABL.kappa

   **type:** Real, optional, default = 0.41
   
   Wall model coefficient.
   
.. input_param:: ABL.surface_roughness_z0

   **type:** Real, optional, default = 0.1
   
   Wall model surface roughness. 
   
.. input_param:: ABL.normal_direction

   **type:** Integer, optional, default = 2
   
   Wall model normal direction. x-direction = 0, y-direction = 1, z-direction = 2.

.. input_param:: ABL.log_law_height

   **type:** Real, optional
   
   Height to evaluate the log law for the wall model. 
   Currently, if this parameter is not specified in the input file, the first half cell height is calculated
   and is used to set the log law height. 
   Therefore the log law height depends on the domain size and number of elements in the normal direction. 
   If this parameter is set to a number the log law is evaluated at a fixed height. 
   Note: currently the fluctuating velocity terms in the shear stress model are only 
   available at the first cell center above the wall. This limitation will be removed soon.
   
   
.. input_param:: ABL.temperature_heights

   **type:** List of Reals, mandatory 
   
   Height(s) in meters at which temperature values are prescribed. 
   
.. input_param:: ABL.temperature_values

   **type:** List of Reals (has to be same length as :input_param:`ABL.temperature_heights`), mandatory
   
   Temperature values in Kelvin at the corresponding :input_param:`ABL.temperature_heights`.
   The temperature below the first height is assumed to be constant and equal to the
   first temperature value. 
   The temperature between values is initialized to have linear variation. 
   The final temperature is constant above the last specified height. 
   
      
.. input_param:: ABL.perturb_velocity

   **type:** Boolean, optional, default = true
   
   If true this flag turns on perturbations to the freestream flow. 
   
.. input_param:: ABL.pertub_ref_height

   **type:** Real, optional, default = 50.0
   
   Reference height for velocity perturbations, 
   perturbations exist below this height and decay above this height.
   Only active when :input_param:`ABL.perturb_velocity` = true.
   
.. input_param:: ABL.Uperiods

   **type:** Real, optional, default = 4.0
   
   Number of sinusoidal waves in x-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.
   
.. input_param:: ABL.Vperiods

   **type:** Real, optional, default = 4.0
   
   Number of sinusoidal waves in y-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.
   
.. input_param:: ABL.deltaU

   **type:** Real, optional, default = 1.0
   
   Amplitude of fluctuations in x-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.
   
.. input_param:: ABL.deltaV

   **type:** Real, optional, default = 1.0
   
   Amplitude of fluctuations in y-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.
   
.. input_param:: ABL.perturb_temperature

   **type:** Boolean, optional, default = false
   
   Perturb temperature field with random fluctuations.
   This feature is currently deactivated. 
	
	
Section: ``Momentum Sources``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: ICNS.source_terms

   **type:** String(s), optional
   
   Activates source terms for the incompressible Navier-Stokes momentum equations. 
   These strings can be entered in any order with a space between each. 
   Current existing source terms for ABL include "BoussinesqBuoyancy", "CoriolisForcing",
   and "ABLForcing". 
   Other source terms include "MMSForcing" and "DensityBuoyancy".
   
.. input_param:: BoussinesqBuoyancy.reference_temperature

   **type:** Real, mandatory
   
   Reference temperature in Kelvin. 
   Values of the temperature field that are less than or greater than this value will 
   cause a buoyancy force in the direction of the gravity vector.
   
.. input_param:: BoussinesqBuoyancy.thermal_expansion_coeff

   **type:** Real, optional, default = ``1/BoussinesqBuoyancy.reference_temperature``
   
   Thermal expansion coefficient, if not specified this value is set to the inverse of the
   :input_param:`BoussinesqBuoyancy.reference_temperature` value.
   
.. input_param:: CoriolisForcing.latitude 

   **type:** Real, mandatory
   
   Latitude in degrees where the Coriolis forcing is computed (northern hemisphere).
   
.. input_param:: CoriolisForcing.rotational_time_period 

   **type:** Real, optional, default = 86400.0
   
   Rotational time period of a day in seconds.
   
.. input_param:: CoriolisForcing.east_vector

   **type:** List of 3 reals, optional, default = 1.0 0.0 0.0
   
   East vector that gives the orientation of the grid w.r.t. to planetary coordinate system.
   This vector is automatically normalized within amr-wind.
   
.. input_param:: CoriolisForcing.north_vector

   **type:** List of 3 reals, optional, default = 0.0 1.0 0.0
   
   North vector that gives the orientation of the grid w.r.t. to planetary coordinate system.
   This vector is automatically normalized within amr-wind.

.. input_param:: ABLForcing.abl_forcing_height 

   **type:** Real, mandatory
   
   height in meters at which the flow is forced to maintain the freestream inflow velocities
   ``incflo.ic_u`` and ``incflo.ic_v``. 

Section: ``Static Refinement``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section is for controlling the static mesh refinement of the grid. 

.. input_param:: tagging.static_refinement 

   **type:** Boolian, optional, default = false
   
   Static refinement with Cartesian-aligned bounding boxes. 
   
.. input_param:: tagging.static_refinement_def

   **type:** String
   
   Static refinement with Cartesian-aligned bounding boxes input file name. 
   
Section: ``Boundary conditions``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
This section controls the boundary conditions. Only non-periodic BC's need to be defined here.

.. input_param:: xlo.type (or ylo.type, zlo.type, xhi.type, yhi.type, zhi.type)

   **type:** String, mandatory if non-periodic
   
   Boundary condition type on the lo (or hi) side of the domain. 
   Current options are: periodic, pressure_inflow, pressure_outflow, mass_inflow, 
   no_slip_wall, slip_wall, and wall_model. 

.. input_param:: xlo.temperature (or ylo.temperature, etc)

   **type:** Real, optional, default = 0.0
   
   Specifies a temperature gradient at the wall, only activated for slip_wall and wall_model BC types. 
   
Section: ``MLMG options``
~~~~~~~~~~~~~~~~~~~~~~~~~~

This section specifies the Multi-Level Multi-Grid (MLMG) options for each type of linear solve. 
There are three types of linear solves performed in amr-wind "diffusion" which is a cell based 
Helmholtz like solve to advance the momentum equations,
"nodal_proj" is a node based pressure projection, and "mac_proj" projects velocities to faces. 
The options are the same for each and the prefix determines which MLMG option is being specified.
Below the diffusion options are described but the same options apply to "nodal_proj" and "mac_proj".
   
.. input_param:: diffusion.mg_verbose

   **type:** Integer, optional, default = 0
   
   Sets the verbosity of the MLMG solver.
   
.. input_param:: diffusion.mg_cg_verbose

   **type:** Integer, optional, default = 0
   
   Sets the verbosity of the bottom solver within MLMG.
   
.. input_param:: diffusion.mg_max_iter

   **type:** Integer, optional, default = 200
   
   Sets the max number of multigrid iterations before MLMG fails and aborts.
  
.. input_param:: diffusion.mg_cg_max_iter

   **type:** Integer, optional, default = 200
   
   Sets the max number of bottom solver iterations.
   
.. input_param:: diffusion.mg_fmg_max_iter 

   **type:** Integer, optional, default = 0
   
   Sets the number of F-cycle MG iterations to perform before switching to V-cycle MG.
   
.. input_param:: diffusion.mg_max_coarsening_level

   **type:** Integer, optional, default = 100
   
   This parameter sets the max number of multigrid coarsening allowed at the lowest amr 
   level to the solver. 
   Typically setting to a large number means coarsen as much as possible until the grid 
   can not be coarsened anymore.
      
.. input_param:: diffusion.mg_max_order

   **type:** Integer, optional, default = 2
   
   Order of the one-sided stencil applied near physical boundaries and fine/coarse boundaries.
   
.. input_param:: diffusion.mg_mg_rtol

   **type:** Real, optional, default = 1.0e-11
   
   Set the relative tolerance for the linear solver
   
.. input_param:: diffusion.mg_mg_atol

   **type:** Real, optional, default = 1.0e-14
   
   Set the absolute tolerance for the linear solver
   
.. input_param:: diffusion.bottom_solver_type

   **type:** String, optional, default = "bicgstab"
   
   Set the bottom solver type. Current bottom solver options 
   include: smoother, bicgstab, cg, bicgcg, cgbicg, hypre, and petsc. 
   The hyper and petsc options will require compiling with those libraries. 
   
