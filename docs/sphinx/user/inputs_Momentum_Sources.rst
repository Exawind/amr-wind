Section: Momentum Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: ICNS.source_terms

   **type:** String(s), optional
   
   Activates source terms for the incompressible Navier-Stokes momentum
   equations. These strings can be entered in any order with a space between
   each. Please consult `AMR-Wind developer documentation
   <https://exawind.github.io/amr-wind/api_docs/group__icns__src.html>`_ for a
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

.. input_param:: ABLForcing.velocity_timetable

   **type:** String, optional
   
   Input file name for table that lists time, wind speed, and wind direction
   of ABL forcing velocity.