.. _amrwind-abl-bndry-io:

AMR-Wind ABL Boundary I/O
=========================

Wind farm simulations typically require ABL inflow conditions. As
such, a precursor ABL simulation is often performed to collect inflow
conditions for the wind farm. AMR-Wind leverages NetCDF to collect ABL
inflow variables at each time step and write the data to a file. This
file can then be read during a wind farm simulation to populate the
inflow.

.. note::

   - Currently, it is only possible to write `xlo` and `ylo` boundaries.

   - Inflow conditions are linearly interpolated between output times.

   - The time for the simulation that is reading the inflow file must be entirely contained within the inflow times.

   - The simulation reading the inflow file must have the same grid resolution at the boundaries.


Generating the inflow file from an ABL simulation
-------------------------------------------------

The following section can be added to the input file to generate an
inflow file during an ABL simulation:

.. code-block:: none

   ABL.bndry_file = "bndry_file.nc"
   ABL.bndry_io_mode = 0
   ABL.bndry_planes = ylo xlo
   ABL.bndry_output_start_time = 2.0
   ABL.bndry_var_names = velocity temperature

In the case of using the OneEqKsgsM84 model the tke field is also needed.

.. code-block:: none

   ABL.bndry_var_names = velocity temperature tke

Using an inflow file in an ABL simulation
-----------------------------------------

The following section can be added to the input file to read an
inflow file to populate the boundary conditions:

.. code-block:: none

   ABL.bndry_file = "../orig/bndry_file.nc"
   ABL.bndry_io_mode = 1
   ABL.bndry_var_names = velocity temperature

Again, In the case of using the OneEqKsgsM84 model the tke field is also needed.

.. code-block:: none

   ABL.bndry_var_names = velocity temperature tke

The boundary conditions need to be adjusted from periodic to inflow/outflow.
The following lines show the changes that need to be made to the input file
for the x coordinate (similar change for y coordinate when needed):

.. code-block:: none

   geometry.is_periodic           = 0 1 0                          # Periodicity x y z (0/1)

   xlo.type = "mass_inflow"
   xlo.density = 1.0
   xhi.type = "pressure_outflow"
  

Inflow file structure
---------------------

The inflow file is written by the NetCDF library in the following structure:

  - Top level contains the data common to all the groups (title,
    dimensions, time)
  - Mid level contains the groups of planes: Each plane (e.g. `xlo`)
    is assigned a group at this level. It contains variables such as
    the plane normal and perpendicular directions.
  - Bottom levels contains groups of AMR levels, i.e. the levels
    associated with each plane: Each AMR level is assigned a group
    (e.g. "level_0", "level_1", etc) containing the variable
    dimensions and the inflow variables (e.g. velocity) associated
    with that level.

For a multi-level file, `ncdump -h <file>` provides:

.. literalinclude:: ./ncdump_bndry_file.txt
   :linenos:


The inflow file can be inspected with Python as such:

.. raw:: html
   :file: ./inspect_abl_io.html
