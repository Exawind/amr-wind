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
   
.. input_param:: ABL.wall_shear_stress_type

   **type:** String, optional, default = "Moeng"

   Wall shear stress model: options include 
   "constant", "local", "Schumann", and "Moeng"

.. input_param:: ABL.three_ComponentForcing

    **type:** Boolean, optional, default = false

    If this flag is true then all three components
    of the coriolis forcing and geostrophic wind
    forcing are included in the ICNS source terms.
    (Default = two horizontal components in the
    coriolis and geostrophic forcing).

.. input_param:: ABL.initial_condition_input_file

   **type:** String, optional, default= ""
    
   File that contains initial conditions for the
   velocity field in netcdf file format.
   This file is expected to have the same dimensions as the simulation.
   Values are passed directly from the file to the velocity field inside the code.
   Only spanwise velocity components are supported. 
