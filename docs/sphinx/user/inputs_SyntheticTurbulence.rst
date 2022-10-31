Section: SyntheticTurbulence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section is for setting turbulence injection parameters.

.. input_param:: SynthTurb.turbulence_file

   **type:** String, required
   
   Name of the netcdf file that contains the data.
   
.. input_param:: SynthTurb.wind_direction

   **type:** Real, optional, default = 270.
   
   The wind direction. 
   
.. input_param:: SynthTurb.grid_location 

   **type:** List of Reals, required
  
   Location of the turbulent grid.

.. input_param:: SynthTurb.mean_wind_type

   **type:** String, required 
  
   The type of profile to use. Options are ConstValue, LinearProfile, and PowerLawProfile.

.. input_param:: SynthTurb.grid_spacing 

   **type:** Real, optional, default = 2.0
  
   The grid spacing used in the box.

.. input_param:: SynthTurb.time_offset

   **type:** Real, optional, default = 0.0
  
   The time offset betweent he data and the simulation.

.. input_param:: SynthTurb.zref

   **type:** Real, required
  
   The default reference height at the center of the turbulence grid

.. input_param:: SynthTurb.shear_exponent

   **type:** Real, required
  
   The shear exponent value.
   
.. input_param:: SynthTurb.uref

   **type:** Real list, required
  
   The reference value of the velocity vector

.. input_param:: SynthTurb.zoffset

   **type:** Real, optional, default = 
  
   The offset in the z direction between the turbulence box and the simulation

.. input_param:: SynthTurb.umin

   **type:** Real, required
  
   The minimum velocity allowed in the turbulence box.

.. input_param:: SynthTurb.umax

   **type:** Real, required
  
   The maximum velocity allowed in the turbulence box
   
