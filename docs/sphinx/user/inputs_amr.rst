Section: amr
~~~~~~~~~~~~~~~~

This section contains input parameters used by the core AMReX mesh data
structure ``AmrCore`` to determine the base mesh and adaptive mesh refinement
strategies. The refinement criteria is specified through :ref:`inputs_tagging`

.. input_param:: amr.n_cell

   **type:** List of 3 integers, mandatory

   The number of cells in the coarset level of AMR hierarchy

.. input_param:: amr.max_level

   **type:** Integer, optional, default: 0

   The maximum AMR level in the refinement hierarchy. Default value is ``0``
   indicating a single mesh level with uniform resolution in the three
   directions.

.. input_param:: amr.max_grid_size

   **type:** Integer, optional, default: 32

   Maximum number of cells at level 0 in each grid.
   There are also options to specify this value in each direction,
   please refer to AMReX documentation.

.. input_param:: amr.blocking_factor

   **type:** Integer, optional, default: 8

   Each grid must be divisible by :input_param:`amr.blocking_factor`.
   There are also options to specify this value in each direction,
   please refer to AMReX documentation.



