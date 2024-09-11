.. _inputs_momentum_sources:
   
Section: Momentum Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: ICNS.source_terms

   **type:** String(s), optional
   
   Activates source terms for the incompressible Navier-Stokes momentum
   equations. These strings can be entered in any order with a space between
   each. Please consult `AMR-Wind developer documentation
   <https://exawind.github.io/amr-wind/api_docs/group__icns__src.html>`_ for a
   comprehensive list of all momentum source terms available. Note that the
   following input arguments specific to each source term will only be active
   if the corresponding source term (the root name) is listed in 
   :input_param:`ICNS.source_terms`.

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
