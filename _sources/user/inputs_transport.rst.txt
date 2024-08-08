Section: transport
~~~~~~~~~~~~~~~~~~~~~~

This section is for setting thermal and momentum diffusivities.

.. input_param:: transport.model

   **type:** String, optional, default = ConstTransport

   Currently this is the only transport model implemented.
   
.. input_param:: transport.viscosity

   **type:** Real, optional, default = 1.0e-5

   Sets the dynamic viscosity
   
.. input_param:: transport.laminar_prandtl 

   **type:** Real, optional, default = 1.0

   Sets the laminar Prandtl number.
   
.. input_param:: transport.turbulent_prandtl 

   **type:** Real, optional, default = 1.0

   Sets the turbulent Prandtl number.
   
