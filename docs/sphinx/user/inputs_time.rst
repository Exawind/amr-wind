.. _inputs_time:

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
   
.. input_param:: time.regrid_start

  **type:** Integer, optional, default = 0; default = start index upon restart

  This user-specified parameter sets the base timestep onwards which the mesh is adaptively
  refined.

.. input_param:: time.plot_start

  **type:** Integer, optional, default = 0; default = start index upon restart

  This user-specified parameter sets the base timestep onwards which the output (plot files)
  are written to the disk.

.. input_param:: time.checkpoint_start

  **type:** Integer, optional, default = 0; default = start index upon restart

  This user-specified parameter sets the base timestep onwards which the checkpoint (restart) 
  files are written to the disk.

.. input_param:: time.use_force_cfl

   **type:** Boolean, optional, default = true

   If this flag is true then the forces (including the pressure gradient) are included
   in the CFL calculation.
   
