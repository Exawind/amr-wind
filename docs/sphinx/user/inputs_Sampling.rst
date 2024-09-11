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

.. input_param:: sampling.output_delay

   **type:** Integer, optional, default = 0

   Specify the output delay (in timesteps) when data sampling and output will begin
   during a simulation. E.g., a delay of 100 will wait until the 100th timestep to 
   check if, according to the output frequency, sampling should be performed and 
   output to disk.

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

.. input_param:: sampling.int_fields

   **type:** List of one or more strings

   List of CFD simulation int fields to sample and output (e.g. mask_cell)

.. input_param:: sampling.derived_fields

   **type:** List of one or more strings

   List of CFD simulation derived fields to sample and output (e.g. mag_vorticity)

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
sampled by specifying the ``offset_vector`` vector along which the the planes are
offset for as many planes as there are entries in the ``offset`` array.

Example::

  sampling.plane1.type          = PlaneSampler
  sampling.plane1.axis1         = 1.0 0.0 0.0
  sampling.plane1.axis2         = 0.0 0.0 1.0
  sampling.plane1.origin        = 0.0 0.0 0.0
  sampling.plane1.num_points    = 10 10
  sampling.plane1.offset_vector = 1.0 0.0 0.0
  sampling.plane1.offsets       = 0.0 2.0 3.0

Illustration of this example:

.. figure:: planesampler.png
   :alt: PlaneSampler
   :width: 800

   Example of sampling on planes.

Sampling at arbitrary locations
````````````````````````````````

The ``ProbeSampler`` allows the user to sample the flow field at arbitrary
locations read from a text file (default: ``probe_locations.txt``).

Example::

  sampling.probe1.type = ProbeSampler
  sampling.probe1.probe_location_file = "probe_locations.txt"

The first line of the file contains the total number of probes for this set.
This is followed by the coordinates (three real numbers), one line per probe.

Sampling on a volume
`````````````````````

The ``VolumeSampler`` samples a 3D volume that starts at ``lo`` and
extends to ``hi``. The resolution in all directions is specified by
``num_points``.

Example::

  sampling.volume1.type        = VolumeSampler
  sampling.volume1.hi        = 3.0 1.0 0.5
  sampling.volume1.lo      = 0.0 0.0 -0.5
  sampling.volume1.num_points  = 30 10 10

Sampling on the air-water interface
```````````````````````````````````

The ``FreeSurfaceSampler`` samples on the air-water interface, and it requires the 
vof (volume-of-fluid) field to be present in order to function. The sample locations
are specified using a grid that starts at ``plane_start`` and
extends to ``plane_end``. The resolution in each direction is specified by
``plane_num_points``. The coordinates of the sampling
locations are determined by the location of the air-water interface in the search
direction, specified by ``search_direction``, and the other coordinates are 
determined by the ``plane_`` parameters. The default search direction parameter
is 2, indicating the samplers will search for the interface along the z-direction. 
Due to this design, it is best to specify a plane that is normal to the intended 
search direction. Another optional parameter is ``num_instances``, which is available
for cases where the interface location is multivalued along the search direction,
such as during wave breaking. This parameter defaults to 1, and the sampler will
automatically select the highest position along the search direction when the interface
location is multivalued.

Example::

  sampling.fs1.type             = FreeSurfaceSampler
  sampling.fs1.plane_start      = 4.0 -1.0 0.0
  sampling.fs1.plane_end        = 0.0 1.0  0.0
  sampling.fs1.plane_num_points = 20 10
