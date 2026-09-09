.. _inputs_turbulence:

Section: turbulence
~~~~~~~~~~~~~~~~~~~

This section is for setting turbulence model parameters

.. input_param:: turbulence.model

   **type:** String, optional, default = Laminar

   Specifies which turbulence model to use, by default "Laminar" is
   chosen (effectively no turbulence model).  Currently the supported
   turbulence models are "Smagorinsky", "AMD", "Kosovic", 
   "OneEqKsgsM84", "KOmegaSST", "KOmegaSSTIDDES" or "KLAxell".

   
.. input_param:: Smagorinsky_coeffs.Cs

   **type:** Real, optional, default = 0.135

   Specifies the coefficient used in the `Smagorinsky` turbulence model. 
   

   

.. input_param:: Kosovic.terrain_model

   **type:** String, optional, default = ``TerrainDrag``

   Terrain physics the ``Kosovic`` model couples to. With terrain the SGS
   viscosity is set to zero inside the terrain and, in the first fluid cells,
   replaced by the log-law value :math:`2 \rho u_*^2 / |\partial U_t/\partial n|`,
   so that the SGS flux across the terrain interface equals the wall stress
   :math:`\rho u_*^2` (the face viscosity is the mean of the cell and the blanked
   neighbor, hence the factor 2).

   - ``TerrainDrag``: uses the ``terrain_blank``, ``terrain_drag``,
     ``terrain_height`` and ``terrainz0`` fields of the
     :ref:`TerrainDrag <inputs_terraindrag>` physics; the log-law value is applied
     in the drag cells with the friction velocity from the cell above at
     :math:`1.5 \Delta z`. Inactive when those fields do not exist.
   - ``ImmersedTerrain``: uses the fields of the
     :ref:`ImmersedTerrain <inputs_immersedterrain>` physics. The viscosity and
     the non-linear term are weighted by one minus the drag weight of the cell
     (:input_param:`ImmersedTerrain.drag_weight`), and in surface cells the
     log-law value is averaged over the same wall patches as the
     ``ImmersedDragForcing`` wall model: the six faces touching a solid
     neighbor, or the surface normal, following
     :input_param:`ImmersedDragForcing.wall_model`,
     :input_param:`ImmersedDragForcing.reference_distance`,
     :input_param:`ImmersedDragForcing.minimum_z0` and
     :input_param:`ImmersedTerrain.solid_threshold`. The tangential speed
     gradient is taken between the surface cell and the wall-side cell across
     the patch, with the same floor of 0.01 1/s as ``TerrainDrag``.
     Requires ``ImmersedTerrain`` in :input_param:`incflo.physics`.
