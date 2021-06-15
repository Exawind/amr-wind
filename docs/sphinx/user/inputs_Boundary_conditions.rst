Section: Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
This section controls the boundary conditions. Only non-periodic BC's need to be defined here.

.. input_param:: xlo.type (or ylo.type, zlo.type, xhi.type, yhi.type, zhi.type)

   **type:** String, mandatory if non-periodic
   
   Boundary condition type on the lo (or hi) side of the domain. 
   Current options are: periodic, pressure_inflow, pressure_outflow, mass_inflow, 
   no_slip_wall, slip_wall, and wall_model. 

.. input_param:: xlo.temperature (or ylo.temperature, etc)

   **type:** Real, optional, default = 0.0
   
   Specifies a temperature gradient at the wall, only activated for slip_wall and wall_model BC types. 
   
