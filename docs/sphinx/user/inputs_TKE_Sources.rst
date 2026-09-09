.. _inputs_tke_sources:

Section: TKE Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: TKE.source_terms

   **type:** String(s), optional
   
These terms are used when the turbulence model includes a transport equation for  turbulent 
kinetic energy. This term has to be set to `KransAxell` for `KLAxell` model.
With terrain, ``KransAxell`` follows :input_param:`KLAxell.terrain_model` to
decide whether the wall and damping forcing is evaluated from the
``TerrainDrag`` or the ``ImmersedTerrain`` fields.
