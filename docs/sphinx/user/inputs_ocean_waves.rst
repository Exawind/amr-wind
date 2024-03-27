Section: Ocean Waves
~~~~~~~~~~~~~~~~~~~~~~~

This section is for setting up wave forcing and relaxation zones.

.. input_param:: OceanWaves.label

   **type:** String

   Name for Ocean Waves instance. As a placeholder, the label 'label' will be used in this documentation.

.. input_param:: OceanWaves.label.type

   **type:** String

   Type of wave to be used. Options include LinearWaves, StokesWaves, and W2AWaves. 
   Linear and stokes waves are monochromatic waves. W2AWaves uses the Waves2AMR library to read 
   and transform HOS-Ocean wave data into the AMR-Wind domain, and this option requires the Waves2AMR 
   library (abbreviated W2A) to be enabled at compile time.

.. input_param:: OceanWaves.label.relax_zone_gen_length

   **type:** Real, optional, default = 4.0

   Length of region, in x-direction, over which the wave profile will be forced. 
   The region begins at the low x boundary and extends in the positive x direction, covering
   the entire width of the domain in the y direction and applying to all liquid at all depths
   (in z). Always on, even if argument is not specified.

.. input_param:: OceanWaves.label.numerical_beach_length

   **type:** Real, optional, default = 8.0

   Length of region, in x-direction, over which the waves will be forced to a flat interface. 
   The region begins at the high x boundary and extends in the negative x direction, covering
   the entire width of the domain in the y direction and applying to all liquid at all depths
   (in z). Always on, unless relax_zone_out_length is included in the input file.

.. input_param:: OceanWaves.label.relax_zone_out_length

   **type:** Real, optional, default = 8.0

   Length of region, in x-direction, over which the wave profile will be forced. This input argument
   replaces, and thereby turns off, the numerical beach.
   The region begins at the high x boundary and extends in the negative x direction, covering
   the entire width of the domain in the y direction and applying to all liquid at all depths
   (in z). Only on when included in the input file.

.. input_param:: OceanWaves.label.zero_sea_level

   **type:** Real, optional, default = 0.0

   Default location (in z) of the liquid surface. Where the water would rest in the absence of waves.

.. input_param:: OceanWaves.label.water_depth

   **type:** Real, optional, default = 0.5

   The depth of the water in the simulation, used for calculating analytical wave profiles (linear 
   and stokes waves). Should be equal to the difference between the zero_sea_level and the low z boundary location.

.. input_param:: OceanWaves.label.timeramp_period

   **type:** Real, optional, default = 2.0

   An initial ramp-up period for the wave forcing. Without specifying a period, the wave 
   forcing will begin at full strength.

.. input_param:: OceanWaves.label.initialize_wave_field

   **type:** Boolean, optional, default = false

   By default, the domain will be initialized with a flat interface; if this option
   is turned on, the wave profile will be initialized over the entire domain. If there is a specified
   relax_zone_out_length, this option is automatically turned on.

The following input arguments are only valid for the LinearWaves and StokesWave wave types:

.. input_param:: OceanWaves.label.wave_length

   **type:** Real, mandatory

   The wave length of the wave profile

.. input_param:: OceanWaves.label.wave_height

   **type:** Real, mandatory

   The amplitude of the wave profile

The following input arguments are only valid for the StokesWave wave type:

.. input_param:: OceanWaves.label.order

   **type:** Integer, mandatory

   The order of the Stokes wave formula being used

The following input arguments are only valid for the W2AWaves wave type:

.. input_param:: OceanWaves.label.HOS_modes_filename

   **type:** String, mandatory

   The name of the modes file, output by HOS-Ocean, in the SWENSE format.

.. input_param:: OceanWaves.label.HOS_init_timestep

   **type:** Integer, optional, default = 0

   The time step in the modes file for the AMR-Wind simulation to start at.

.. input_param:: OceanWaves.label.HOS_init_time

   **type:** Integer, optional, default = 0

   The physical time in the modes file for the AMR-Wind simulation to start at.
   This argument is only active if HOS_init_timestep is omitted. AMR-Wind will pick the
   time step in the modes closest to the specified time.

.. input_param:: OceanWaves.label.number_interp_points_in_z

   **type:** Integer, mandatory

   When Waves2AMR converts mode data to spatial data, the z locations to perform the transformation
   must be chosen. AMR-Wind does this by creating a geometric series from the water surface
   to the lower z boundary, starting with a small spacing between points near the surface and 
   then expanding downward. This argument is the total number of z locations where wave modes will be
   transformed to inform the wave velocity field. After the transformation step, the data is then interpolated
   to the AMR-Wind mesh. More points means better resolution over the whole depth, but more points
   also means more work for the solver.

.. input_param:: OceanWaves.label.interp_spacing_at_surface

   **type:** Real, mandatory

   The physical spacing between interpolation points at the water surface.
   This is the most influential parameter for the geometric series that dictates the location of the points in z.
   This should be set to near the mesh spacing in z around the water surface.

.. input_param:: OceanWaves.label.number_interp_above_surface

   **type:** Integer, optional, default = 1

   The number of points placed above the mean water surface for the velocity transformation process. The spacing
   between the points above the surface is equal to the interp_spacing_at_surface. When setting this value, the wave height
   should be considered so that velocity can be accurately computed for portions of the waves above the mean surface.