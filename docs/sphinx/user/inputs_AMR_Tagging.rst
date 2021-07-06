.. _inputs_tagging:

Section: AMR Tagging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section manages the various mesh refinement criteria that can be used to
activate either static or adaptive mesh refinement during simulations. The
parameters are read from the prefix ``tagging`` and can contain different types
of tagging logic. Note that this section is only active if
:input_param:`amr.max_level` is greater than zero. Regridding interval is controlled by :input_param:`time.regrid_interval` .

Example::

  tagging.labels = s1 f1 g1
  tagging.s1.type = CartBoxRefinement
  tagging.s1.static_refinement_def = static_box.txt

  tagging.f1.type = FieldRefinement
  tagging.f1.field_name = density
  tagging.f1.grad_error = 0.1. 0.1 0.1

  tagging.g1.type = GeometryRefinement
  tagging.g1.shapes = c1 b1

  tagging.g1.c1.type = cylinder
  tagging.g1.c1.start = 500.0 500.0 250.0
  tagging.g1.c1.end = 500.0 500.0 750.0
  tagging.g1.c1.outer_radius = 300.0
  tagging.g1.c1.inner_radius = 275.0

  tagging.g1.b1.type = box
  tagging.g1.b1.origin = 300.0 150.0 250.0
  tagging.g1.b1.xaxis =  450.0 600.0 0.0
  tagging.g1.b1.yaxis =  -150.0 100.0 0.0
  tagging.g1.b1.zaxis = 0.0 0.0 500.0

Each section must contain the keyword ``type`` that is one of the refinement types:

======================== ===================================================================
``CartBoxRefinement``    Nested refinement using Cartesian boxes
``FieldRefinement``      Refinement based on error metric for field or its gradient
``OversetRefinement``    Refinement around fringe/field interface
``GeometryRefinement``   Refinement using geometric shapes
``QCriterionRefinement`` Refinement using Q-Criterion
``VorticityMagRefinement``  Refinement using vorticity
======================== ===================================================================

.. input_param:: tagging.labels

   **type:** List of one or more names

   Labels indicate a list of prefixes for different types of refinement criteria
   active during the simulation.

The parameters for the subsections are determined by the type of refinement being performed.

Refinement using Cartesian boxes
````````````````````````````````

``CartBoxRefinement`` allows refining boxes (aligned with the principal axes).

Example::

   tagging.labels = static
   tagging.static.type = CartBoxRefinement
   tagging.static.static_refinement_def = static_box.txt

.. input_param:: tagging.CartBoxRefinement.static_refinement_def

   **type:** String, required

   The text file that contains a list of bounding boxes used to perform
   refinement at various levels.

Refinement using field error criteria
`````````````````````````````````````

Example::

  tagging.f1.type = FieldRefinement
  tagging.f1.field_name = density
  tagging.f1.grad_error = 0.1. 0.1 0.1

.. input_param:: tagging.FieldRefinement.field_name

   **type:** String, required

   The name of the field used to tag cells

.. input_param:: tagging.FieldRefinement.field_error

   **type:** Vector<Real>, optional

   List of field error values at each level. The user must specify a value for
   each level desired.

.. input_param:: tagging.FieldRefinement.grad_error

   **type:** Vector<Real>, optional

   List of gradient error values at each level. The user must specify a value for
   each level desired.

Refinement using geometry
`````````````````````````

This section controls refinement using pre-defined geometric shapes. Currently,
two options are supported: 1. ``box`` -- refines the region inside a hexahedral
block, and 2. ``cylinder`` -- refines the region inside a cylindrical block.

.. input_param:: tagging.GeometryRefinement.shapes

   **type:** List of strings, required

   Names of the input subsections that define specific geometries for refinement.

.. input_param:: tagging.GeoemtryRefinement.level

   **type:**  Integer, optional, default: -1

   If ``level`` is provided and is greater than or equal to 0, then the
   refinement based on geometries defined for this section is only performed at
   that level.

.. input_param:: tagging.GeometryRefinement.min_level

   **type:**  Integer, optional, default: 0

   If ``level`` is not specified, then this option specifies the minimum level
   where this refinement is active.

.. input_param:: tagging.GeometryRefinement.max_level

   **type:**  Integer, optional, default: ``mesh.maxLevel()``

   If ``level`` is not specified, then this option specifies the maximum level
   where this refinement is active.

Note that the specification of ``level`` overrides, ``min_level`` and
``max_level`` specifications. This can be used to control the different levels
where refinement regions are active.

Example::

  tagging.g1.type = GeometryRefinement
  tagging.g1.shapes = b1 b2
  tagging.g1.level = 0
  tagging.g1.b1.type = box
  tagging.g1.b1.origin = 300.0 150.0 250.0
  tagging.g1.b1.xaxis =  450.0 600.0 0.0
  tagging.g1.b1.yaxis =  -150.0 100.0 0.0
  tagging.g1.b1.zaxis = 0.0 0.0 500.0
  tagging.g1.b2.type = box
  tagging.g1.b2.origin = 600.0 350.0 250.0
  tagging.g1.b2.xaxis =  50.0 30.0 0.0
  tagging.g1.b2.yaxis =  -50.0 60.0 0.0
  tagging.g1.b2.zaxis = 0.0 0.0 500.0

  tagging.g2.type = GeometryRefinement
  tagging.g2.shapes = c1
  tagging.g2.level = 1
  tagging.g2.c1.type = cylinder
  tagging.g2.c1.start = 500.0 500.0 250.0
  tagging.g2.c1.end = 500.0 500.0 750.0
  tagging.g2.c1.outer_radius = 300.0
  tagging.g2.c1.inner_radius = 275.0


This example defines two different refinement definitions acting on level 0 and
1 respectively. The refinement at level 0 (``g1``) contains two box regions,
whereas the refinement at level 1 (``g2``) only contains one cylinder
definition.

**Refinement using hexahedral block definitions**

To perform ``box`` refinement, the user specifies the ``origin`` of the box and
three vectors: ``xaxis, yaxis, zaxis`` that defines the directions and the
extents of the hexahedral block. Denoting :math:`\mathbf{O}` as origin vector
and :math:`\mathbf{x}`, :math:`\mathbf{y}` and :math:`\mathbf{z}` as the three
vectors given by the user, the position vectors of the eight corners of the
hexahedral box are given by

.. math::

   \mathbf{x}_0 &= \mathbf{O} && \mathbf{x}_4 &= \mathbf{O} + \mathbf{z} \\
   \mathbf{x}_1 &= \mathbf{O} + \mathbf{x} && \mathbf{x}_5 &= \mathbf{O} + \mathbf{z} + \mathbf{x} \\
   \mathbf{x}_2 &= \mathbf{O} + \mathbf{x} + \mathbf{y} \qquad && \mathbf{x}_6 &= \mathbf{O} + \mathbf{z} + \mathbf{x} + \mathbf{y} \\
   \mathbf{x}_3 &= \mathbf{O} + \mathbf{y} && \mathbf{x}_7 &= \mathbf{O} + \mathbf{z} + \mathbf{y} \\



**Refinement using cylindrical block definitions**

The axis and the extents along the axis are defined by two position vectors
``start`` and ``end``. The radial extent is specified by ``outer_radius``. An
optional ``inner_radius`` can be specified to restrict tagging to an annulus
between the inner and outer radii.

Refinement using Q-Criterion
`````````````````````````````````````

Example::

  tagging.qc1.type = QCriterionRefinement
  tagging.qc1.nondim = false
  tagging.qc1.values = 10.0 20.0 20.0

.. input_param:: tagging.QCriterionRefinement.nondim

   **type:** Boolean, optional, default = true

   Boolean determining if the dimensional or non-dimensional form 
   of Q-criterion should be used. Dimensional version may require 
   modifying values depending on physical scales. For the non-dimensional 
   form positive thresholds indicate regions where the rotational strength is 
   larger than the shear rate strength. A threshold of unity indicates 
   that the rotational strength is equal to the background shear strength. 
   
.. input_param:: tagging.QCriterionRefinement.values

   **type:** Vector<Real>, optional

   List of Q-criterion values at each level.
   If the absolute value of Q-criterion exceeds this value
   the cell is tagged for refinement.
   The user must specify a value for each level desired.

