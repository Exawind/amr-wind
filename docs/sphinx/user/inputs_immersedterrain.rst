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
- ``terrain_surface`` (Real, 3 components): terrain height at the cell center
  and the slopes :math:`\partial h/\partial x` and :math:`\partial h/\partial y`.
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
   center. ``distance_function`` uses a smooth hyperbolic tangent of the signed distance from
   the cell center to the surface, with the width set by
   :input_param:`ImmersedTerrain.smoothing_length`.

.. input_param:: ImmersedTerrain.smoothing_length

   **type:** Real, optional, default = 1.0

   Smoothing length for ``distance_function`` blanking in units of the vertical
   cell size.

.. input_param:: ImmersedTerrain.solid_threshold

   **type:** Real, optional, default = 0.5

   A neighboring cell is treated as a wall when its terrain fraction is at or
   above this value. Used to classify surface cells and, in
   ``ImmersedDragForcing``, to decide which faces of a surface cell carry the
   wall model.

.. input_param:: ImmersedTerrain.implicit_projection

   **type:** Boolean, optional, default = false

   Apply the immersed drag implicitly through the nodal and MAC projections
   instead of as an explicit source term. The terrain then behaves as a fluid of
   density :math:`\rho (1 + \beta C \Delta t)` in the pressure solve, with
   :math:`C = C_d / \Delta z` and :math:`C_d` taken from
   :input_param:`ImmersedDragForcing.drag_coefficient`, so that the pressure
   gradient produces no velocity inside the terrain. Without it the projection
   re-injects :math:`\Delta t \nabla p / \rho` inside the body every step and the
   velocity residual inside the terrain decreases only linearly with the time
   step. When active, ``ImmersedDragForcing`` skips its explicit drag term and
   applies only the wall model. Declares the field ``terrain_drag_rate``.

.. input_param:: ImmersedTerrain.interface_diffusion

   **type:** String, optional, default = ``none``

   Treatment of the diffusive flux across the terrain interface. Without it the
   diffusion operator computes a flux :math:`\mu_\mathrm{eff} (u_k - 0)/\Delta_f`
   at every fluid/solid face, because the interior is at rest, which is a wall stress
   with the wrong length scale and, with a turbulence model, double counts the stress
   supplied by the wall model.

   - ``none``: current behavior.
   - ``block``: the face coefficient is multiplied by :math:`\min(1-\beta_L, 1-\beta_R)`,
     which removes the molecular and SGS flux across faces touching the terrain so that the
     wall model alone carries the wall stress and heat flux. Intended for the turbulent pathway.
   - ``no_slip``: on faces between a fluid cell and a cell with fraction at or above
     :input_param:`ImmersedTerrain.solid_threshold`, the coefficient is scaled by
     :math:`\Delta_f / d_1` so that the discrete flux equals :math:`\mu u_k / d_1`, the flux to
     a no-slip wall at the true position. On the bottom face :math:`d_1 = z_k - h` from the
     terrain height (clamped to :math:`[0.1, 1] \Delta z`); on other faces the wall is taken at
     the face, :math:`d_1 = \Delta_f/2`. Intended for laminar flow with a finite viscosity.

   The factors are stored in the face fields ``terrain_diffusion_xf/yf/zf`` and applied
   to every equation that uses the shared diffusion operator (momentum, temperature, TKE,
   passive scalars).

.. input_param:: ImmersedTerrain.drag_weight

   **type:** String, optional, default = ``fraction``

   Weight of the immersed drag in partially filled cells, also used as the complement of
   the wall-model weight. ``fraction`` uses the terrain fraction :math:`\beta` (drag
   :math:`\beta`, wall model :math:`1-\beta`). ``center`` treats a cell whose center is
   inside the terrain (:math:`\beta \ge` :input_param:`ImmersedTerrain.solid_threshold`)
   as fully solid and any other partial cell as a fluid cell that receives the full wall
   model at its true distance and no drag. With the ``terrain_height`` or
   ``surface_normal`` wall models this is the cut-cell configuration whose wall-model error
   is first order per cell; with ``fraction`` a mostly fluid partial cell is also damped
   by the drag, which over-damps the first layer on slopes by a factor independent of
   the mesh. Applies to the momentum and temperature sources and to the implicit
   projection rate.
