.. _inputs_immersedterrain:

Section: ImmersedTerrain
~~~~~~~~~~~~~~~~~~~~~~~~

These parameters are active when ``ImmersedTerrain`` is included in
:input_param:`incflo.physics`. ImmersedTerrain is the successor to
:ref:`TerrainDrag <inputs_terraindrag>`: instead of a binary blanking it
stores the fraction of each cell occupied by terrain, so slopes are not
represented as a staircase. It is used together with the
``ImmersedDragForcing`` momentum source (see
:ref:`inputs_momentum_sources`).

ImmersedTerrain declares the following fields, each with one ghost cell:

- ``terrain_fraction`` (Real): fraction of the cell volume occupied by terrain,
  0 in fluid, 1 inside terrain, in between for cells cut by the surface.
- ``terrain_mask`` (Int): 0 for fluid, 1 for cells fully inside the terrain,
  2 for surface cells, i.e. partially filled cells and fluid cells that share
  a face with a mostly solid cell on any of the six sides. The wall model acts
  in surface cells. The mask can be used directly with
  ``FieldRefinement``: ``field_error = 1.5`` tags only the surface band and
  ``field_error = 0.5`` tags surface and solid cells.
- ``terrain_surface`` (Real, 3 components): terrain height at the cell centre
  and the slopes dh/dx and dh/dy.
- ``terrain_roughness`` (Real): aerodynamic roughness length.

.. input_param:: ImmersedTerrain.terrain_file

   **type:** String, optional, default = ``terrain.amrwind``

   Input file for terrain height data, in the same flat-grid format as
   :input_param:`TerrainDrag.terrain_file`.

.. input_param:: ImmersedTerrain.roughness_file

   **type:** String, optional, default = ``terrain.roughness``

   Input file for roughness-length data. If this file is missing or cannot be
   opened, the roughness is set from
   :input_param:`ImmersedTerrain.uniform_roughness`.

.. input_param:: ImmersedTerrain.uniform_roughness

   **type:** Real, optional, default = 0.1

   Uniform roughness length used when the roughness file is not available.

.. input_param:: ImmersedTerrain.blanking_method

   **type:** String, optional, default = ``volume_fraction``

   How the terrain fraction of a cell is computed. ``volume_fraction`` uses the
   fraction of the cell column that lies below the terrain height at the cell
   centre. ``distance_function`` uses a smooth tanh of the signed distance from
   the cell centre to the surface, with the width set by
   :input_param:`ImmersedTerrain.smoothing_length`.

.. input_param:: ImmersedTerrain.smoothing_length

   **type:** Real, optional, default = 1.0

   Smoothing length for ``distance_function`` blanking in units of the vertical
   cell size.

.. input_param:: ImmersedTerrain.solid_threshold

   **type:** Real, optional, default = 0.5

   A neighbouring cell is treated as a wall when its terrain fraction is at or
   above this value. Used to classify surface cells and, in
   ``ImmersedDragForcing``, to decide which faces of a surface cell carry the
   wall model.
