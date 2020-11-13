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
``Boundary conditions`` Boundary condition types and gradients
``MLMG options``        Multi-Level Multi-Grid Linear solver options
``Tagging``             Static and dynamic refinement options
``Sampling``            Data probes to sample field data during simulations
======================= ============================================================

This section documents the parameters available within each section.

.. note::

   - Boolean flags (``true/false``) can also be indicated using integers in the text file
     and uses the convention ``0 = False`` and ``1 = True``.

   - Quotes around strings are optional.

   - If an input is repeated only the last one is used.

Section: geometry
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

Section: amr
~~~~~~~~~~~~~~~~

This section contains input parameters used by the core AMReX mesh data
structure ``AmrCore`` to determine the base mesh and adaptive mesh refinement
strategies. The refinement criteria is specified through :ref:`inputs_tagging`

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

   Each grid must be divisible by :input_param:`amr.blocking_factor`.
   There are also options to specify this value in each direction,
   please refer to AMReX documentation.

Section: time
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
   
Section: io
~~~~~~~~~~~~~~~~~

This section deals with parameters that control input/output to the simulation.

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
   
   
Section: incflo
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

.. input_param:: incflo.post_processing

   **type:** List of strings, optional

   When present, this parameter contains list of sections to be read with
   specific post-postprocessing actions. Currently, the code supports
   :ref:`Sampling <inputs_sampling>`.

   ::

     incflo.post_processing     = sampling
     sampling.output_frequency  = 5
     sampling.labels            = line1 line2
     sampling.fields            = velocity
     sampling/line1.type        = LineSampler
     sampling/line1.num_points  = 21
     sampling/line1.start       = 250.0 250.0 10.0
     sampling/line1.end         = 250.0 250.0 210.0
     sampling/line2.type        = LineSampler
     sampling/line2.num_points  = 21
     sampling/line2.start       = 500.0 500.0 10.0
     sampling/line2.end         = 500.0 500.0 210.0

   In the above example, the code will read the parameters with keyword
   ``sampling`` to initialize user-defined probes.
   
Section: transport
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
   
Section: turbulence
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
   
Section: ABL
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

.. input_param:: ABL.theta_amplitude

   **type:** Real, optional, default = 0.8 K

   Amplitude of the temperature perturbations added to the initial field. Only
   active when :input_param:`ABL.perturb_temperature` is true.

.. input_param:: ABL.cutoff_height

   **type:** Real, optional, default = domain height

   Height below which temperature perturbations are added

.. input_param:: ABL.random_gauss_mean

   **type:** Real, optional, default = 0.0

   Mean for the Gaussian random number generator

.. input_param:: ABL.random_gauss_var

   **type:** Real, optional, default = 1.0

   Variance for the Gaussian random number generator

	
.. input_param:: ABL.bndry_file

   **type:** String, optional, default = ""

   NetCDF-4 file name for ABL inflow

.. input_param:: ABL.bndry_io_mode

   **type:** Int, optional, default = -1

   IO mode (0=output, 1=input)

.. input_param:: ABL.bndry_planes

   **type:** String, optional, default = ""

   IO planes for ABL inflow

.. input_param:: ABL.bndry_output_start_time

   **type:** Real, optional, default = 0.0

   Time at which to start ABL inflow output

.. input_param:: ABL.bndry_var_names

   **type:** String, optional, default = ""

   Variables for IO for ABL inflow


Section: Momentum Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: ICNS.source_terms

   **type:** String(s), optional
   
   Activates source terms for the incompressible Navier-Stokes momentum
   equations. These strings can be entered in any order with a space between
   each. Please consult `AMR-Wind developer documentation
   <https://exawind.github.io/amr-wind/group__icns__src.html>`_ for a
   comprehensive list of all momentum source terms available.

.. input_param:: BoussinesqBuoyancy.reference_temperature

   **type:** Real, mandatory
   
   Reference temperature :math:`\theta_\mathrm{ref}` in Kelvin.
   Values of the temperature field that are less than or greater than this value will 
   cause a buoyancy force in the direction of the gravity vector.
   
.. input_param:: BoussinesqBuoyancy.thermal_expansion_coeff

   **type:** Real, optional, default :math:`\beta = 1 / \theta_\mathrm{ref}`
   
   Thermal expansion coefficient, if not specified this value is set to the inverse of the
   :input_param:`BoussinesqBuoyancy.reference_temperature` value.
   
.. input_param:: CoriolisForcing.latitude 

   **type:** Real, mandatory
   
   Latitude in degrees where the Coriolis forcing is computed. Positive values
   indicate northern hemisphere.
   
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

.. input_param:: GeostrophicForcing.geostrophic_wind

   **type:** List of 3 reals, optional

   The user has to choose between GeostrophicForcing and ABLForcing. 
   CoriolisForcing input must be present when using GeostrophicForcing.
   These checks are not enforced for now.
   
.. input_param:: ABLForcing.abl_forcing_height

   **type:** Real, mandatory
   
   Height in meters at which the flow is forced to maintain the freestream
   inflow velocities specified through :input_param:`incflo.velocity`.

Section: Static Refinement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section is for controlling the static mesh refinement of the grid. 

.. input_param:: tagging.static_refinement 

   **type:** Boolean, optional, default = false
   
   Static refinement with Cartesian-aligned bounding boxes. 
   
.. input_param:: tagging.static_refinement_def

   **type:** String
   
   Static refinement with Cartesian-aligned bounding boxes input file name. 
   
Section: Boundary conditions
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
   
Section: MLMG options
~~~~~~~~~~~~~~~~~~~~~~~~~~

This section specifies the Multi-Level Multi-Grid (MLMG) options for each type
of linear solve. There are three types of linear solves performed in amr-wind
"diffusion" which is a cell based Helmholtz like solve to advance the momentum
equations, "nodal_proj" is a node based pressure projection, and "mac_proj"
projects velocities to faces. The options are the same for each and the prefix
determines which MLMG option is being specified. Below the diffusion options are
described but the same options apply to "nodal_proj" and "mac_proj". It is also
possible to specify diffusion solver options for specific equations such as
temperature, to do that use "temperature_diffusion" as your prefix.

**Linear operator options**
   
.. input_param:: diffusion.max_coarsening_level

   **type:** Integer, optional, default = 100
   
   This parameter sets the max number of multigrid coarsening allowed at the lowest amr 
   level to the solver. 
   Typically setting to a large number means coarsen as much as possible until the grid 
   can not be coarsened anymore.
      
.. input_param:: diffusion.max_order

   **type:** Integer, optional, default = 2
   
   Order of the one-sided stencil applied near physical boundaries and fine/coarse boundaries.

**MLMG options**

.. input_param:: diffusion.verbose

   **type:** Integer, optional, default = 0

   Sets the verbosity of the MLMG solver.

.. input_param:: diffusion.maxiter

   **type:** Integer, optional, default = 200

   Sets the max number of multigrid iterations. If :input_param:`do_fixed_iters`
   is set to True, then AMReX will not abort if specified tolerance is not met
   after max iterations, otherwise it will abort.

.. input_param:: diffusion.do_fixed_iters

   **type:** Boolean, optional, default = true

   If ``true``, then AMReX will not abort if the specified tolerance is not met
   even after :input_param:`diffusion.maxiter` iterations have completed.

.. input_param:: diffusion.mg_rtol

   **type:** Real, optional, default = 1.0e-11
   
   Set the relative tolerance for the linear solver
   
.. input_param:: diffusion.mg_atol

   **type:** Real, optional, default = 1.0e-14
   
   Set the absolute tolerance for the linear solver

.. input_param:: diffusion.fmg_maxiter

   **type:** Integer, optional, default = 0

   Sets the number of F-cycle MG iterations to perform before switching to V-cycle MG.

.. input_param:: diffusion.num_pre_smooth

   **type:** Integer, optional, default = 2

   Number of pre smoothing steps

.. input_param:: diffusion.num_post_smooth

   **type:** Integer, optional, default = 2

   Number of post smoothing steps

.. input_param:: diffusion.num_final_smooth

   **type:** Integer, optional, default = 8

   Number of final smoother steps applied

.. input_param:: diffusion.num_bottom_smooth

   **type:** Integer, optional, default = 0

   Number of smoother steps applied during bottom solve.

**Bottom solver options**
   
.. input_param:: diffusion.bottom_solver

   **type:** String, optional, default = "bicgstab"
   
   Set the bottom solver type. Current bottom solver options 
   include: smoother, bicgstab, cg, bicgcg, cgbicg, hypre, and petsc. 
   The hyper and petsc options will require compiling with those libraries.

.. input_param:: diffusion.bottom_verbose

   **type:** Integer, optional, default = 0

   Sets the verbosity of the bottom solver within MLMG.

.. input_param:: diffusion.bottom_rtol

   **type:** Real, optional, default = 1.0e-4

   Set the relative tolerance for the bottom solver for convergence.

.. input_param:: diffusion.bottom_atol

   **type:** Real, optional, default = -1.0

   Set the absolute tolerance for the bottom solve. Setting a negative number
   disables absolute tolerance check.

.. input_param:: diffusion.bottom_maxiter

   **type:** Integer, optional, default = 200

   Maximum number of iterations for the bottom solver

.. input_param:: diffusion.hypre_interface

   **type:** String, optional, default = ``ij``

   The hypre interface to use when :input_param:`diffusion.bottom_solver` is set
   to ``hypre``. Valid choices are: ``ij``, ``semi_structured``, and
   ``structured``.

.. input_param:: diffusion.hypre_namespace

   The ParmParse ``prefix`` where the hypre options must be read from for this
   solver. For example, to set hypre options for NodalProjector

   ..
      nodal_proj.hypre_namespace = "nodal_proj.hypre"
      nodal_proj.hypre.hypre_solver = GMRES
      nodal_proj.hypre.hypre_preconditioner = BoomerAMG

.. _inputs_tagging:

Section: AMR Tagging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section manages the various mesh refinement criteria that can be used to
activate either static or adaptive mesh refinement during simulations. The
parameters are read from the prefix ``tagging`` and can contain different types
of tagging logic. Note that this section is only active if
:input_param:`amr.max_level` is greater than zero. Regridding interval is controlled by :input_param:`time.regrid_interval` .

Example::

  tagging.labels = s1 f1 g1
  tagging/s1.type = CartBoxRefinement
  tagging/s1.static_refinement_def = static_box.txt

  tagging/f1.type = FieldRefinement
  tagging/f1.field_name = density
  tagging/f1.grad_error = 0.1. 0.1 0.1

  tagging/g1.type = GeometryRefinement
  tagging/g1.shapes = c1 b1

  tagging/g1/c1.type = cylinder
  tagging/g1/c1.start = 500.0 500.0 250.0
  tagging/g1/c1.end = 500.0 500.0 750.0
  tagging/g1/c1.outer_radius = 300.0
  tagging/g1/c1.inner_radius = 275.0

  tagging/g1/b1.type = box
  tagging/g1/b1.origin = 300.0 150.0 250.0
  tagging/g1/b1.xaxis =  450.0 600.0 0.0
  tagging/g1/b1.yaxis =  -150.0 100.0 0.0
  tagging/g1/b1.zaxis = 0.0 0.0 500.0

Each section must contain the keyword ``type`` that is one of the refinement types:

======================== ===================================================================
``CartBoxRefinement``    Nested refinement using Cartesian boxes
``FieldRefinement``      Refinement based on error metric for field or its gradient
``OversetRefinement``    Refinement around fringe/field interface
``GeometryRefinement``   Refinement using geometric shapes
``QCriterionRefinement`` Refinement using Q-Criterion
``VorticityRefinement``  Refinement using vorticity
======================== ===================================================================

.. input_param:: tagging.labels

   **type:** List of one or more names

   Labels indicate a list of prefixes for different types of refinement criteria
   active during the simulation.

The parameters for the subsections are determined by the type of refinement being performed.

Refinement using Cartesian boxes
````````````````````````````````

``CartBoxRefinement`` allows refining boxes (aligned with the principal axes).

Example::

   tagging.labels = static
   tagging/static.type = CartBoxRefinement
   tagging/static.static_refinement_def = static_box.txt

.. input_param:: tagging.CartBoxRefinement.static_refinement_def

   **type:** String, required

   The text file that contains a list of bounding boxes used to perform
   refinement at various levels.

Refinement using field error criteria
`````````````````````````````````````

Example::

  tagging/f1.type = FieldRefinement
  tagging/f1.field_name = density
  tagging/f1.grad_error = 0.1. 0.1 0.1

.. input_param:: tagging.FieldRefinement.field_name

   **type:** String, required

   The name of the field used to tag cells

.. input_param:: tagging.FieldRefinement.field_error

   **type:** Vector<Real>, optional

   List of field error values at each level. The user must specify a value for
   each level desired.

.. input_param:: tagging.FieldRefinement.grad_error

   **type:** Vector<Real>, optional

   List of gradient error values at each level. The user must specify a value for
   each level desired.

Refinement using geometry
`````````````````````````

This section controls refinement using pre-defined geometric shapes. Currently,
two options are supported: 1. ``box`` -- refines the region inside a hexahedral
block, and 2. ``cylinder`` -- refines the region inside a cylindrical block.

.. input_param:: tagging.GeometryRefinement.shapes

   **type:** List of strings, required

   Names of the input subsections that define specific geometries for refinement.

.. input_param:: tagging.GeoemtryRefinement.level

   **type:**  Integer, optional, default: -1

   If ``level`` is provided and is greater than or equal to 0, then the
   refinement based on geometries defined for this section is only performed at
   that level.

.. input_param:: tagging.GeometryRefinement.min_level

   **type:**  Integer, optional, default: 0

   If ``level`` is not specified, then this option specifies the minimum level
   where this refinement is active.

.. input_param:: tagging.GeometryRefinement.max_level

   **type:**  Integer, optional, default: ``mesh.maxLevel()``

   If ``level`` is not specified, then this option specifies the maximum level
   where this refinement is active.

Note that the specification of ``level`` overrides, ``min_level`` and
``max_level`` specifications. This can be used to control the different levels
where refinement regions are active.

Example::

  tagging/g1.type = GeometryRefinement
  tagging/g1.shapes = b1 b2
  tagging/g1.level = 0
  tagging/g1/b1.type = box
  tagging/g1/b1.origin = 300.0 150.0 250.0
  tagging/g1/b1.xaxis =  450.0 600.0 0.0
  tagging/g1/b1.yaxis =  -150.0 100.0 0.0
  tagging/g1/b1.zaxis = 0.0 0.0 500.0
  tagging/g1/b2.type = box
  tagging/g1/b2.origin = 600.0 350.0 250.0
  tagging/g1/b2.xaxis =  50.0 30.0 0.0
  tagging/g1/b2.yaxis =  -50.0 60.0 0.0
  tagging/g1/b2.zaxis = 0.0 0.0 500.0

  tagging/g2.type = GeometryRefinement
  tagging/g2.shapes = c1
  tagging/g2.level = 1
  tagging/g2/c1.type = cylinder
  tagging/g2/c1.start = 500.0 500.0 250.0
  tagging/g2/c1.end = 500.0 500.0 750.0
  tagging/g2/c1.outer_radius = 300.0
  tagging/g2/c1.inner_radius = 275.0


This example defines two different refinement definitions acting on level 0 and
1 respectively. The refinement at level 0 (``g1``) contains two box regions,
whereas the refinement at level 1 (``g2``) only contains one cylinder
definition.

**Refinement using hexahedral block definitions**

To perform ``box`` refinement, the user specifies the ``origin`` of the box and
three vectors: ``xaxis, yaxis, zaxis`` that defines the directions and the
extents of the hexahedral block. Denoting :math:`\mathbf{O}` as origin vector
and :math:`\mathbf{x}`, :math:`\mathbf{y}` and :math:`\mathbf{z}` as the three
vectors given by the user, the position vectors of the eight corners of the
hexahedral box are given by

.. math::

   \mathbf{x}_0 &= \mathbf{O} && \mathbf{x}_4 &= \mathbf{O} + \mathbf{z} \\
   \mathbf{x}_1 &= \mathbf{O} + \mathbf{x} && \mathbf{x}_5 &= \mathbf{O} + \mathbf{z} + \mathbf{x} \\
   \mathbf{x}_2 &= \mathbf{O} + \mathbf{x} + \mathbf{y} \qquad && \mathbf{x}_6 &= \mathbf{O} + \mathbf{z} + \mathbf{x} + \mathbf{y} \\
   \mathbf{x}_3 &= \mathbf{O} + \mathbf{y} && \mathbf{x}_7 &= \mathbf{O} + \mathbf{z} + \mathbf{y} \\



**Refinement using cylindrical block definitions**

The axis and the extents along the axis are defined by two position vectors
``start`` and ``end``. The radial extent is specified by ``outer_radius``. An
optional ``inner_radius`` can be specified to restrict tagging to an annulus
between the inner and outer radii.

Refinement using Q-Criterion
`````````````````````````````````````

Example::

  tagging/qc1.type = QCriterionRefinement
  tagging/qc1.nondim = false
  tagging/qc1.values = 10.0 20.0 20.0

.. input_param:: tagging.QCriterionRefinement.nondim

   **type:** Boolean, optional, default = true

   Boolean determining if the dimensional or non-dimensional form 
   of Q-criterion should be used. Dimensional version may require 
   modifying values depending on physical scales. For the non-dimensional 
   form positive thresholds indicate regions where the rotational strength is 
   larger than the shear rate strength. A threshold of unity indicates 
   that the rotational strength is equal to the background shear strength. 
   
.. input_param:: tagging.QCriterionRefinement.values

   **type:** Vector<Real>, optional

   List of Q-criterion values at each level.
   If the absolute value of Q-criterion exceeds this value
   the cell is tagged for refinement.
   The user must specify a value for each level desired.

.. _inputs_sampling:
   
Section: Sampling
~~~~~~~~~~~~~~~~~~~~~

This section controls data-sampling (post-processing) actions supported within
AMR-wind. Note that while the input parameters use the keyword ``sampling``, the
actual keyword is determined by the labels provided to
:input_param:`incflo.post_processing`. So for example, if
``incflo.post_processing = my_sampling``, then the options must be prefixed with
``my_sampling.``.

.. input_param:: sampling.output_frequency

   **type:** Integer, optional, default = 100

   Specify the output frequency (in timesteps) when data sampling is performed
   and output to disk.

.. input_param:: sampling.output_format

   **type:** String, optional, default = "native"

   Specify the format of the data outputs. Currently the code supports the
   following formats

   ``native``
       AMReX particle binary format

   ``ascii``
       AMReX particle ASCII format. Note, this can have significant impact
       on performance and must be used for debugging only.

.. input_param:: sampling.labels

   **type:** List of one or more names

   Labels indicate the names of the different types of samplers (e.g., line,
   plane, probes) that are used to sample data from the flow field.

   For example, if the user uses

   Example::

      sampling.labels = line1 plane1 probe1

   Then the code expects to read ``sampling/line1, sampling/plane1,
   sampling/probe1`` sections to determine the specific sampling probe information.

.. input_param:: sampling.fields

   **type:** List of one or more strings

   List of CFD simulation fields to sample and output

The individual sampling types are documented below

Sampling along a line
``````````````````````

The ``LineSampler`` allows the user to sample the flow-field along a line
defined by ``start`` and ``end`` coordinates with ``num_points`` equidistant
nodes.

Example::

  sampling/line1.type       = LineSampler
  sampling/line1.num_points = 21
  sampling/line1.start      = 250.0 250.0 10.0
  sampling/line1.end        = 250.0 250.0 210.0

Sampling on one or more planes
```````````````````````````````

The ``PlaneSampler`` samples the flow-field on two-dimensional planes defined by
two axes: ``axis1`` and ``axis2`` with the bottom corner located at ``origin``
and is divided into equally spaced nodes defined by the two entries in
``num_points`` vector. Multiple planes parallel to the reference planes can be
sampled by specifying the ``normal`` vector along which the the planes are
offset for as many planes as there are entries in the ``offset`` array.

Example::

  sampling/plane1.type        = PlaneSampler
  sampling/plane1.axis1       = 0.0 1.0 0.0
  sampling/plane1.axis2       = 0.0 0.0 1.0
  sampling/plane1.origin      = 0.0 0.0 0.0
  sampling/plane1.num_points  = 10 10
  sampling/plane1.normal      = 1.0 0.0 0.0
  sampling/plane1.offsets     = -10.0 0.0 10.0

Sampling at arbitrary locations
````````````````````````````````

The ``ProbeSampler`` allows the user to sample the flow field at arbitrary
locations read from a text file (default: ``probe_locations.txt``).

Example::

  sampling/probe1.type = ProbeSampler
  sampling/probe1.probe_location_file = "probe_locations.txt"

The first line of the file contains the total number of probes for this set.
This is followed by the coordinates (three real numbers), one line per probe.
