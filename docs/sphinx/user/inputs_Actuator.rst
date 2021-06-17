
.. _inputs_actuator:

Section: Actuator
~~~~~~~~~~~~~~~~~~

This section controls  the actuator type models. This includes the actuator 
disk and line models. The prefix is the label set in 
``incflo.physics``. For example
``incflo.physics = FreeStream Actuator``
Actuator models are meant to simulate aerodynamic objects by using body forces 
in the momentum equation. 
There are capabilities to simulate fixed wings as actuator lines and wind
turbines as actuator disks and actuator line models.


.. input_param:: Actuator.labels

   **type:** String, mandatory
   
   This string is used as an identifier for the current actuator.


.. input_param:: Actuator.type

   **type:** String, mandatory
   
   This string identifies the type of actuator to use. The ones currently
   supported are: ``TurbineFastLine``, ``TurbineFastDisk``, and 
   ``FixedWingLine``.

FixedWingLine
"""""""""""""

Example for ``FixedWingLine``::

   incflo.physics = FreeStream Actuator 
   Actuator.labels = F1 
   Actuator.type = FixedWingLine 
   Actuator.FixedWingLine.num_points = 21 
   Actuator.FixedWingLine.epsilon = 3.0 3.0 3.0 
   Actuator.FixedWingLine.pitch = 4.0 
   Actuator.FixedWingLine.span_locs = 0.0 1.0 
   Actuator.FixedWingLine.chord = 2.0 2.0 
   Actuator.FixedWingLine.airfoil_table = DU21_A17.txt 
   Actuator.FixedWingLine.airfoil_type = openfast 
   Actuator.F1.start = 0.0 -4.0 0.0 
   Actuator.F1.end = 0.0 4.0 0.0 
   Actuator.F1.output_frequency = 10 
   ICNS.source_terms = ActuatorForcing 

.. input_param:: Actuator.FixedWingLine.num_points

   **type:** int, mandatory
   
   This is the number of actuator points along the wing to be used in the 
   simulation.

.. input_param:: Actuator.FixedWingLine.epsilon

   **type:** List of 3 real numbers, mandatory
   
   This is the value of epsilon in the chord, thicness and spanwise directions.

.. input_param:: Actuator.FixedWingLine.epsilon_chord

   **type:** List of 3 real numbers, optional
   
   This is the value of epsilon/chord. This value will be used to compute 
   epsilon as a function of the chord at every actuator point. A value of 
   epsilon / chord ~ 0.2 is recommended for an optimal representation of the 
   blade aerodynamics. When this variable is specified, the code will choose
   the maximum value between ``epsilon_chord * chord`` and ``epsilon`` for
   every actuator point.

.. input_param:: Actuator.FixedWingLine.pitch

   **type:** Real number, optional
   
   This is the pitch angle of the blade in degrees. All coordinates will be 
   pitched by this angle. In the case of a fixed wing, this would be the angle
   of attack of the wing with respect to the inflow velocity.

.. input_param:: Actuator.FixedWingLine.span_locs

   **type:** List of real numbers, mandatory

   These are non-dimensional span locations from 0 to 1. These locations are
   used to specify the chord values at avery span location of the blade.

.. input_param:: Actuator.FixedWingLine.chord

   **type:** List of real numbers, mandatory
   
   These are the chord values at every span location. The length of this array 
   needs to be the same length as ``span_locs``.

.. input_param:: Actuator.FixedWingLine.airfoil_table

   **type:** String, mandatory
   
   This is the name of the file that contains the lookup table for lift and drag
   coefficients.

.. input_param:: Actuator.FixedWingLine.airfoil_type

   **type:** String, mandatory
   
   This is the type of airfoil table lookup. The currently supported options are
   ``openfast`` and ``text``.

.. input_param:: Actuator.F1.start

   **type:** List of 3 real numbers, mandatory
   This is the starting point of the wing where the first actuator point will be.

.. input_param:: Actuator.F1.end

   **type:** List of 3 real numbers, mandatory
   This is the end point of the wing where the last actuator point will be.

.. input_param:: Actuator.F1.output_frequency

   **type:** int, optional
   
   This is how often to write actuator output. The default is ``10``.


TurbineFastLine
"""""""""""""""

Example for ``TurbineFastLine``::

   incflo.physics = FreeStream Actuator
   Actuator.labels = WTG01
   Actuator.type = TurbineFastLine
   Actuator.TurbineFastLine.rotor_diameter = 126.0
   Actuator.TurbineFastLine.hub_height = 90.0
   Actuator.TurbineFastLine.num_points_blade = 64
   Actuator.TurbineFastLine.num_points_tower = 12
   Actuator.TurbineFastLine.epsilon = 10.0 10.0 10.0
   Actuator.TurbineFastLine.epsilon_tower = 5.0 5.0 5.0
   Actuator.TurbineFastLine.openfast_start_time = 0.0
   Actuator.TurbineFastLine.openfast_stop_time = 1.0
   Actuator.TurbineFastLine.nacelle_drag_coeff = 0.0
   Actuator.TurbineFastLine.nacelle_area = 0.0
   Actuator.TurbineFastLine.output_frequency = 10
   Actuator.TurbineFastLine.density = 1.225
   Actuator.WTG01.base_position = 5.0191 0. -89.56256
   Actuator.WTG01.openfast_input_file = "fast_inp/nrel5mw.fst"
   ICNS.source_terms = ActuatorForcing

.. input_param:: Actuator.TurbineFastLine.rotor_diameter

   **type:** Real number, required
   
   This is the rotor diameter of the turbine to be simulated.

.. input_param:: Actuator.TurbineFastLine.hub_height

   **type:** Real number, required
   
   This is the hub height of the turbine.

.. input_param:: Actuator.TurbineFastLine.num_points_blade

   **type:** int, required
   
   This the number of actuator points along the blades.

.. input_param:: Actuator.TurbineFastLine.num_points_tower

   **type:** int, required
   
   This is the number of actuator points along the tower.

.. input_param:: Actuator.TurbineFastLine.openfast_start_time

   **type:** Real, required
   
   This is the time at which to start the openfast simulation.

.. input_param:: Actuator.TurbineFastLine.openfast_stop_time

   **type:** Real, required
   
   This is the time at which to stop the openfast run.

.. input_param:: Actuator.TurbineFastLine.nacelle_drag_coeff 

   **type:** Real, optional
   
   This is the drag coefficient of the nacelle. If this and the area of the 
   nacelle are specified, a value of epsilon for the nacelle is computed that
   would provide an optimal momentum thickness of the wake.
   This is also used to correct the sampled velocity at the location of the
   nacelle actuator point.

.. input_param:: Actuator.TurbineFastLine.nacelle_area

   **type:** Real, optional, default=0
   
   This is the frontal area of the nacelle which is used to compute the force.

.. input_param:: Actuator.TurbineFastLine.output_frequency

   **type:** int, optional, default=10
   
   This is how often to write actuator output. 

.. input_param:: Actuator.TurbineFastLine.density

   **type:** Real, optional
   
   This is the density of the fluid specified in openfast. This is used to 
   non-dimensionalize the forces from openfast.

.. input_param:: Actuator.WTG01.openfast_input_file

   **type:** String, required
   
   This is the name of the openfast input file with all the turbine information.






