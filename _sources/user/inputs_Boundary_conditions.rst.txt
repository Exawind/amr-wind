Section: Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
This section controls the boundary conditions. Only non-periodic BC's need to be defined here.

.. input_param:: xlo.type (or ylo.type, zlo.type, xhi.type, yhi.type, zhi.type)

   **type:** String, mandatory if non-periodic
   
   Boundary condition type on the lo (or hi) side of the domain. 
   Current options are: periodic, pressure_inflow, pressure_outflow, mass_inflow, 
   no_slip_wall, slip_wall, symmetric_wall and wall_model.

.. input_param:: xlo.temperature (or ylo.temperature, etc)

   **type:** Real, optional, default = 0.0
   
   Specifies a temperature gradient at the wall, only activated for slip_wall and wall_model BC types. 

Mass inflow boundary conditions
```````````````````````````````

For the mass_inflow boundary condition type, the user can specify constant values in the input file::

  xlo.type = "mass_inflow"
  xlo.velocity = 1.0 1.0 0.0
  xlo.temperature = 300.0

In addition to being able to specify constant values, the user has access to a variety of included complex boundary conditions (e.g., LinearProfile, PowerLawProfile, Rankine) that are activated in the input file. For example::

  xlo.type = "mass_inflow"
  xlo.density = 1.0
  xlo.velocity = 5.0 5.0 0.0
  xlo.velocity.inflow_type = LinearProfile

If the user wants to define their own boundary conditions, this is done by editing `CustomScalar` and `CustomVelocity` source and header files in the `udfs` folder. `CustomScalar` is used for scalar fields and `CustomVelocity` is used for velocity fields. These can then be activated in the input file as such::

  xlo.type = "mass_inflow"
  xlo.temperature.inflow_type = CustomScalar
  CustomScalar.foo = 1.0
  xlo.velocity.inflow_type = CustomVelocity
  CustomVelocity.foo = 1.0

They do not both need to be defined at the same time. It is the user's responsibility to ensure that the source files are appropriately edited for their use case. Examples of how these files can be edited are found through comparison of the other mass_inflow functions in the `udfs` folder.
