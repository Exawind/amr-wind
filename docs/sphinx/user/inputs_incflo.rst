Section: incflo
~~~~~~~~~~~~~~~~~~~

This section deals with parameters that mostly determine how amr-wind is run such 
as initial conditions and discretization options.

.. input_param:: incflo.physics

   **type:** List of strings

   Specify a string or a list of strings for each type of physics to initialize and simulate.
   Physics is additive and more than one type of physics may be used.
   Current implemented physics are FreeStream, ABL, Actuator, RayleighTaylor, BoussinesqBubble, 
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
   
.. input_param:: incflo.godunov_type

   **type:** String, optional, default = ppm

   Specifies which Godunov scheme to use. Options include ``plm``, ``ppm``, 
   ``ppm_nolim``, ``weno_js``, and ``weno_z``
   
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
   
