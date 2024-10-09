.. _inputs_incflo:

Section: incflo
~~~~~~~~~~~~~~~

This section deals with parameters that mostly determine how amr-wind is run such 
as initial conditions and discretization options.

.. input_param:: incflo.physics

   **type:** List of strings

   Specify a string or a list of strings for each type of physics to initialize and simulate.
   Physics is additive and more than one type of physics may be used.
   Current implemented physics include FreeStream, SyntheticTurbulence, ABL, Actuator, RayleighTaylor, BoussinesqBubble, TaylorGreenVortex, and ScalarAdvection (which is an example of using a passive scalar advection).
   For multiphase simulations, the MultiPhase physics must be specified, and for forcing wave profiles into the domain, the OceanWaves physics must be specified as well.
   For immersed boundary forcing method TerrainDrag must be specified and the folder should include a `terrain.amrwind` file. 
   
.. input_param:: incflo.density

   **type:** Real

   Specify reference density. 
   For the most part if :input_param:`incflo.constant_density` = true then `incflo.density` sets a constant density everywhere.
   Refer to the field initializer for your chosen :input_param:`incflo.physics` for how `incflo.density` is used.
   
.. input_param:: incflo.velocity

   **type:** List of three Real numbers, sometimes mandatory depending on the :input_param:`incflo.physics`

   Specify reference velocity in the "x-", "y-", "z-direction". 
   Refer to the field initializer for your chosen :input_param:`incflo.physics` for how `incflo.velocity` is used.
   In the context of ABL flows, this argument specifies the initial bulk velocity as well as the target velocity
   for ABL forcing terms, unless reference values come from a file instead. Many field initializers do not use this
   input argument.
   
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

   **type:** Boolean, optional, default = true

   Specifies which advection scheme to use: either Godunov (true) or method of lines (false). 
   Godunov the default approach and has many advantages over the method of lines (MOL): better accuracy,
   stability at larger CFL numbers, and greater computational efficiency. Setting this argument to false is 
   not recommended, and active use or development relying on the method of lines (MOL) is very sparse.
   
.. _inputs_incflo_advection:

.. input_param:: incflo.godunov_type

   **type:** String, optional, default = weno_z

   Specifies which Godunov scheme to use. Options include ``plm``, ``ppm``, 
   ``ppm_nolim``, ``weno_js``, and ``weno_z``
   
.. input_param:: incflo.godunov_use_forces_in_trans

   **type:** Boolean, optional, default = false

   Specifies if body forces are included in the transverse velocity prediction.
   Note: only used when :input_param:`incflo.use_godunov` = true.
   
.. _inputs_incflo_diffusion:

.. input_param:: incflo.diffusion_type

   **type:** Integer, optional, default = 2

   Determines how the diffusion term is handled when updating the momentum equations. 
   A value of 0 is explicit diffusion and all diffusion terms are moved to the right hand side 
   (warning this carries with it a more stringent CFL restriction), 
   a value of 1 is Crank-Nicolson and diffusion terms are on both the left and right hand sides,
   and a value of 2 (default) is a fully implicit diffusion where the entire diffusion term is handled on the left hand side.
   
.. input_param:: incflo.post_processing

   **type:** List of strings, optional

   When present, this parameter contains list of sections to be read with
   specific post-postprocessing actions. Currently, the code supports
   :ref:`Sampling <inputs_sampling>`, :ref:`KineticEnergy <inputs_ke>`,
   :ref:`Enstrophy <inputs_enst>` and :ref:`Averaging <inputs_averaging>`

   ::

     incflo.post_processing     = sampling ke enst
     sampling.type              = Sampling
     sampling.output_frequency  = 5
     sampling.labels            = line1 line2
     sampling.fields            = velocity
     sampling.line1.type        = LineSampler
     sampling.line1.num_points  = 21
     sampling.line1.start       = 250.0 250.0 10.0
     sampling.line1.end         = 250.0 250.0 210.0
     sampling.line2.type        = LineSampler
     sampling.line2.num_points  = 21
     sampling.line2.start       = 500.0 500.0 10.0
     sampling.line2.end         = 500.0 500.0 210.0
     ke.type                    = KineticEnergy
     ke.output_frequency        = 2

   In the above example, the code will read the parameters with keyword
   ``sampling`` to initialize user-defined probes.
   
