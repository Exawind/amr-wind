.. _inputs_sampling:
   
Section: Sampling
~~~~~~~~~~~~~~~~~

This section controls data-sampling (post-processing) actions supported within
AMR-wind. Note that while the input parameters use the keyword ``sampling``, the
actual keyword is determined by the labels provided to
:input_param:`incflo.post_processing`. So for example, if
``incflo.post_processing = my_sampling``, then the options must be prefixed with
``my_sampling.``.

.. input_param:: sampling.output_frequency

   **type:** Integer, optional, default = 100

   Specify the output frequency (in timesteps) when data sampling is performed
   and output to disk.

.. input_param:: sampling.output_format

   **type:** String, optional, default = "native"

   Specify the format of the data outputs. Currently the code supports the
   following formats

   ``native``
       AMReX particle binary format

   ``ascii``
       AMReX particle ASCII format. Note, this can have significant impact
       on performance and must be used for debugging only.
       
   ``netcdf``
       This is the preferred output format and requires linking to the
       netcdf library. If netcdf is linked to AMR-Wind and output format 
       is not specified then netcdf is chosen by default.

.. input_param:: sampling.labels

   **type:** List of one or more names

   Labels indicate the names of the different types of samplers (e.g., line,
   plane, probes) that are used to sample data from the flow field.

   For example, if the user uses

   Example::

      sampling.labels = line1 lidar1 plane1 probe1

   Then the code expects to read ``sampling.line1, sampling.plane1,
   sampling.probe1`` sections to determine the specific sampling probe information.

.. input_param:: sampling.fields

   **type:** List of one or more strings

   List of CFD simulation fields to sample and output

The individual sampling types are documented below

Sampling along a line
``````````````````````

The ``LineSampler`` allows the user to sample the flow-field along a line
defined by ``start`` and ``end`` coordinates with ``num_points`` equidistant
nodes.

Example::

  sampling.line1.type       = LineSampler
  sampling.line1.num_points = 21
  sampling.line1.start      = 250.0 250.0 10.0
  sampling.line1.end        = 250.0 250.0 210.0

Sampling along a line moving in time (virtual lidar)
``````````````````````````````````````````````````````

The ``LidarSampler`` allows the user to sample the flow-field along a line
defined by ``origin`` and spanning to ``length`` 
with ``num_points`` equidistant nodes.
Location of the line is given by the time histories 
``azimuth_table`` and ``elevation_table``.
Angles are given in degrees with 0 azimuth and 0 elevation being the 
x direction. Lidar measurements may also be collected at a constant location
by specifying only one entry to the tables.

Example::

  sampling.lidar1.type            = LidarSampler
  sampling.lidar1.num_points      = 21
  sampling.lidar1.origin          = 250.0 250.0 10.0
  sampling.lidar1.length          = 500.0
  sampling.lidar1.time_table      = 0 10.0
  sampling.lidar1.azimuth_table   = 0 90.0
  sampling.lidar1.elevation_table = 0 45.0

Sampling on one or more planes
```````````````````````````````

The ``PlaneSampler`` samples the flow-field on two-dimensional planes defined by
two axes: ``axis1`` and ``axis2`` with the bottom corner located at ``origin``
and is divided into equally spaced nodes defined by the two entries in
``num_points`` vector. Multiple planes parallel to the reference planes can be
sampled by specifying the ``normal`` vector along which the the planes are
offset for as many planes as there are entries in the ``offset`` array.

Example::

  sampling.plane1.type        = PlaneSampler
  sampling.plane1.axis1       = 0.0 1.0 0.0
  sampling.plane1.axis2       = 0.0 0.0 1.0
  sampling.plane1.origin      = 0.0 0.0 0.0
  sampling.plane1.num_points  = 10 10
  sampling.plane1.normal      = 1.0 0.0 0.0
  sampling.plane1.offsets     = -10.0 0.0 10.0

Sampling at arbitrary locations
````````````````````````````````

The ``ProbeSampler`` allows the user to sample the flow field at arbitrary
locations read from a text file (default: ``probe_locations.txt``).

Example::

  sampling.probe1.type = ProbeSampler
  sampling.probe1.probe_location_file = "probe_locations.txt"

The first line of the file contains the total number of probes for this set.
This is followed by the coordinates (three real numbers), one line per probe.


