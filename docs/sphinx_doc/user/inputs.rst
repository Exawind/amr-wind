.. _amrwind-inputs:

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
   :emphasize-lines: 6, 11, 15, 25, 37, 45


Input file reference
---------------------

The AMR-Wind input file is organized in the following sections

==================== ============================================================
Section              Description
==================== ============================================================
``geometry``         Computational domain information
``amr``              Mesh refinement controls
``time``             Simulation time controls
``incflo``           CFD algorithm and physics controls
``abl``              Atmospheric boundary layer (ABL) controls
==================== ============================================================

This section documents the parameters available within each section.

.. note::

   Boolean flags (``true/false``) are indicated using integers in the text file
   and uses the convention ``0 = False`` and ``1 = True``.

Section: ``geometry``
~~~~~~~~~~~~~~~~~~~~~

This section deals with inputs related to the problem domain.

.. input_param:: geometry.prob_lo

   **type:** List of 3 real numbers, mandatory

   The coordiates of *lower corner* of the computational domain bounding box.

.. input_param:: geometry.prob_hi

   **type:** List of 3 real numbers, mandatory

   The coordinates of the *upper corner* of the computational domain bounding box.

.. input_param:: geometry.is_periodic

   **type:** List of 3 integers, mandatory

   Flags indicating whether the flow is periodic in the ``x``, ``y``, or ``z``
   directions respectively.

Section: ``amr``
~~~~~~~~~~~~~~~~

This section mostly contains inputs the adaptive mesh refinement capabilities of
:program:`amr_wind`.

.. input_param:: amr.n_cell

   **type:** List of 3 integers, mandatory

   The number of cells in the coarset level of AMR hierarchy

.. input_param:: amr.max_level

   **type:** Integer, optional, default: 0

   The maximum AMR level in the refinement hierarchy. Default value is ``0``
   indicating a single mesh level with uniform resolution in the three
   directions.

Section: ``time``
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

   **type:** integer

   The number of timesteps to run before termination. See also
   :input_param:`time.stop_time`.

.. input_param:: time.fixed_dt

   **type:** Real number

   Fixed timestep size (in seconds) used to advance the simulation. If this
   parameter is negative, then :input_param:`time.cfl` is used to determine the
   adaptive timestep during the simulation.

.. input_param:: time.cfl

   **type:** Real number

   The maximum CFL allowed during the simulation. When a positive value is
   provided, adaptive timestepping algorithm is used to ensure that the maximum
   CFL condition is not violated while using the largest allowable timestep to
   advance the simulation. If the user specifies a negative value, then
   :input_param:`time.fixed_dt` is used to advance the simulation with fixed
   timestep.

.. input_param:: time.init_shrink

   **type:** Real number, optional, default = 0.1

   ``init_shrink`` indicates a reduction factor applied to the timestep size
   during initialization phased used to perform the initial iterations to
   determine pressure. This factor must be a value greater than zero but less
   than or equal to 1.0.

.. input_param:: time.regrid_interval

   **type:** integer

   If :input_param:`amr.max_level` is greater than zero, this parameter
   indicates the frequency (in timesteps) at which the mesh is adaptively
   refined based on various user-specified criteria. If this value is negative,
   the mesh is only refined once during initialization and remains constant for
   the rest of the simulation.

.. input_param:: time.plot_interval

   **type:** integer

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which outputs (plot files) are written to disk.

.. input_param:: time.checkpoint_interval

   **type:** integer

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which checkpoint (restart) files are written to disk.
