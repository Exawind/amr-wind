.. _inputs_momentum_sources:
   
Section: Momentum Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: ICNS.source_terms

   **type:** String(s), optional
   
   Activates source terms for the incompressible Navier-Stokes momentum
   equations. These strings can be entered in any order with a space between
   each. Please consult the :doc:`../doxygen/html/index` for a
   comprehensive list of all momentum source terms available. Note that the
   following input arguments specific to each source term will only be active
   if the corresponding source term (the root name) is listed in 
   :input_param:`ICNS.source_terms`.

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
   This vector is automatically normalized within AMR-Wind.
   
.. input_param:: CoriolisForcing.north_vector

   **type:** List of 3 reals, optional, default = 0.0 1.0 0.0
   
   North vector that gives the orientation of the grid w.r.t. to planetary coordinate system.
   This vector is automatically normalized within AMR-Wind.

.. input_param:: GeostrophicForcing.geostrophic_wind

   **type:** List of 3 reals, optional

   The user has to choose between GeostrophicForcing and ABLForcing. 
   CoriolisForcing input must be present when using GeostrophicForcing.
   These checks are not enforced for now.

.. input_param:: GeostrophicForcing.geostrophic_wind_timetable

   **type:** String, optional
   
   Input file name for table that lists time in seconds, wind speed 
   in meters per second, and horizontal wind direction in degrees of the Geostrophic 
   forcing velocity. Each line in the file should be a sequence of 
   three floats specifying the inputs in that order (e.g., 0.0 8.0 -5.0). If this
   argument is present, the :input_param:`GeostrophicForcing.geostrophic_wind`
   will be ignored. Note that the code expects there to be a single-line header
   at the beginning of the geostrophic wind timetable file; if no header exists, 
   the first line of data will be ignored.
   
.. input_param:: ABLForcing.abl_forcing_height

   **type:** Real, mandatory
   
   Height in meters at which the flow is forced to maintain the freestream
   inflow velocities specified through :input_param:`incflo.velocity`.

.. input_param:: ABLForcing.velocity_timetable

   **type:** String, optional
   
   Input file name for table that lists time in seconds, wind speed 
   in meters per second, and horizontal wind direction in degrees of the ABL 
   forcing velocity. Each line in the file should be a sequence of 
   three floats specifying the inputs in that order (e.g., 0.0 8.0 -5.0). If this
   argument is present, the :input_param:`incflo.velocity` argument
   will be ignored. Note that the code expects there to be a single-line header
   at the beginning of the velocity timetable file; if no header exists, the first 
   line of data will be ignored.

.. input_param:: ABLForcing.forcing_timetable_output_file

   **type:** String, optional
   
   Output file name for writing the ABL forcing vector to a text file over the course
   of a simulation. This output is primarily intended for replicating the ABL forcing
   from a precursor simulation in a subsequent inflow-outflow simulation by providing 
   an input file for Body Forcing. The output file will contain the time and three vector
   components of the force.

.. input_param:: ABLForcing.forcing_timetable_frequency

   **type:** Int, optional
   
   The interval of timesteps for writing to the forcing timetable output file. The default
   is 1, i.e., writing every step, which is also the default of the boundary plane output feature.

.. input_param:: ABLForcing.forcing_timetable_start_time

   **type:** Real, optional
   
   The start time for writing to the forcing timetable output file. The default is 0.

.. input_param:: ABLForcing.abl_forcing_off_height

   **type:** Real, required for multiphase simulations with ABL
   
   This parameter indicates the vertical distance above the water level that the ABL
   forcing term should be turned off. This tuning parameter is used to avoid applying 
   the ABL forcing to ocean waves. This is not used when the volume fraction field (vof)
   is not present in the simulation.

.. input_param:: ABLForcing.abl_forcing_ramp_height

   **type:** Real, required for multiphase simulations with ABL
   
   This parameter indicates the vertical distance above the water level and the "off height"
   that the ABL forcing term should ramp up from zero to full strength. This is not used
   when the volume fraction field (vof) is not present in the simulation.

.. input_param:: ABLForcing.abl_forcing_band

   **type:** Real, optional for multiphase simulations with ABL
   
   This parameter is an additional safeguard against applying ABL forcing within the waves.
   This specifies the number of computational cells in a band around the air-water interface
   that the ABL forcing should be deactivated. While the other arguments relate to the height coordinate
   within the domain, this argument is relative to the actual position of water in the simulation.
   The default value is 2.

.. input_param:: BodyForce.type

   **type:** String, optional
   
   The type of body force being used. The default is uniform_constant, which applies a single constant
   force vector over the entire domain. Other available types are height_varying, oscillatory, and
   uniform_constant.

.. input_param:: BodyForce.magnitude

   **type:** List of 3 reals, conditionally mandatory
   
   The force vector to be applied as a body force. This argument is mandatory for uniform_constant 
   (default) and oscillatory body force types.

.. input_param:: BodyForce.angular_frequency

   **type:** Real, conditionally mandatory
   
   The angular frequency to be used for applying sinusoidal time variation to the body force. This 
   argument is mandatory for the oscillatory body force type and is only active for the oscillatory type.

.. input_param:: BodyForce.bodyforce_file

   **type:** String, conditionally mandatory
   
   The text file for specifying the body force vector as a function of height. This text file must contain
   heights (z coordinate values), force components in x, and force components in y. This argument is mandatory for
   the height_varying body force type and is only active for the height_varying type.

.. input_param:: BodyForce.uniform_timetable_file

   **type:** String, conditionally mandatory
   
   The text file for specifying the body force vector as a function of time. This text file must contain
   times, force components in x, force components in y, and force components in z. This argument is mandatory for
   the uniform_timetable body force type and is only active for the uniform_timetable type.  Note that the code 
   expects there to be a single-line header at the beginning of the uniform timetable file; if no header exists, 
   the first line of data will be ignored.

.. input_param:: DragForcing.drag_coefficient

   **type:** Real, optional

   This value specifies the coefficient for the forcing term in the immersed boundary forcing method. It is currently
   recommended to use the default value to avoid initial numerical stability.

.. input_param:: DragForcing.is_laminar

   **type:** Boolean, optional, default = false

   This term is intended to indicate whether the simulation is laminar or not. For laminar
   cases, the drag coefficient is reduced.

.. input_param:: DragForcing.wave_model_inviscid_form_drag

   **type:** Boolean, optional, default = false

   This input file option turns on or off an inviscid model for the form drag of waves in the domain. 
   The formulation of this model is adapted from the Moving Surface Drag (MOSD) model developed by
   `Ayala et al (2024) <https://doi.org/10.1007/s10546-024-00884-8>`_.
   
   When the OceanWaves physics module is active, and the volume fraction variable ("vof") is not in the simulation,
   DragForcing will represent ocean waves as moving terrain. This is automatic and independent of DragForcing
   input arguments. When waves are represented as moving terrain and
   there is sufficient mesh resolution to resolve the shape of the wave, the blanking of cells performed
   by the DragForcing routine will naturally introduce the form drag of the waves into the flow. However,
   when the waves are not sufficiently resolved, such as when the wave amplitude is less than the cell height,
   the analytical model for the form drag, activated by setting this option to true, can be used to compensate
   for the lack of resolution. Therefore, this option should remain set to false except in scenarios
   when the form drag is known to be under-resolved.


The following arguments are influential when ``GravityForcing`` is included in :input_param:`ICNS.source_terms`.

   .. input_param:: ICNS.use_perturb_pressure

   **type:** Boolean, optional, default = false
   
   When this option is off, the GravityForcing term is simply :math:`g`, which becomes
   :math:`\rho g` when included in the momentum equation. By activating this option,
   the momentum term applied by GravityForcing will become :math:`(\rho - \rho_0) g`,
   where :math:`rho_0` is some constant reference density profile. The reference density field
   can be created by either MultiPhase physics or anelastic ABL physics. By using the
   reference density, the pressure field seen by the solver is represented as a
   perturbation from a reference pressure field, enabling pressure_outflow boundary
   conditions to better handle certain flows, e.g., those with equilibrium pressure gradients
   parallel to the outflow plane.

   .. input_param:: ICNS.reconstruct_true_pressure

   **type:** Boolean, optional, default = false
   
   This option is only valid when the perturbational pressure form is being used, i.e.,
   :input_param:`ICNS.use_perturb_pressure` = true. Reconstructing the true pressure
   adds back the reference pressure profile to obtain the full pressure after the
   pressure solve has been performed. This makes no difference to the flow evolution,
   but it changes the field available for post-processing or coupling to overset solvers.

Rayleigh Damping is a momentum source term that forces the flow near the boundaries.
All of the default arguments deal with the upper boundary (in z), damping the flow
toward a target velocity. There are also additional arguments to use similar forcing
in the lateral directions.

.. input_param:: RayleighDamping.time_scale

   **type:** Real, required

   Time scale of the damping. The smaller the time scale (and closer to the time step size),
   the stronger and more immediate the damping will be.

.. input_param:: RayleighDamping.length_complete_damping

   **type:** Real, required

   Beginning at the top boundary, the complete damping region uses damping at full strength, i.e. "complete damping".
   This length designates how far the complete damping region extends from the top boundary.

.. input_param:: RayleighDamping.length_sloped_damping

   **type:** Real, required

   Beginning at the bottom of the complete damping region, the sloped damping region gradually
   transitions from complete damping to no damping at all. This length designates how far the
   sloped damping region extends from the bottom of the complete damping region.

.. input_param:: RayleighDamping.reference_velocity

   **type:** Vector<Real>, required

   The target velocity that the damping forces the flow toward.

.. input_param:: RayleighDamping.force_coord_direction

   **type:** Vector<Int>, optional, default = [1, 1, 1]

   This parameter determines which coordinate directions of velocity will be forced.
   The default is to force every direction, but often only forcing in the z direction
   is desired, which can be achieved by setting this parameter to 0 0 1.

.. input_param:: RayleighDamping.lateral_damping_start

   **type:** Real, optional

   There is optional damping in the lateral directions as well, which can be turned on
   with this parameter and the ones that follow. This parameter sets the coordinate in z
   above which the lateral damping (in x and y) is active. By default, this is set to a
   very high value in order to turn it off when not set. Lateral damping only applies
   to the z component of the velocity.

.. input_param:: RayleighDamping.vertical_cutoff

   **type:** Real, optional

   The upper bound (in z) of the region where lateral damping is applied. This argument
   is required when :input_param:`RayleighDamping.lateral_damping_start` is included.

.. input_param:: RayleighDamping.west_damping_length

   **type:** Real, optional, default = 0

   Length of region starting from the west (-x) boundary to apply lateral damping term.

.. input_param:: RayleighDamping.east_damping_length

   **type:** Real, optional, default = 0

   Length of region starting from the east (+x) boundary to apply lateral damping term.

.. input_param:: RayleighDamping.south_damping_length

   **type:** Real, optional, default = 0

   Length of region starting from the south (-y) boundary to apply lateral damping term.

.. input_param:: RayleighDamping.north_damping_length

   **type:** Real, optional, default = 0

   Length of region starting from the north (+y) boundary to apply lateral damping term.
