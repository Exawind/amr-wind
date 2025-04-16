.. _inputs_multiphase:

Section: MultiPhase
~~~~~~~~~~~~~~~~~~~~

This section is for the MultiPhase physics module, which is required for flows that solve air and water together using the volume-of-fluid method.

.. input_param:: MultiPhase.density_fluid1

   **type:** Real, optional, default = 10.

   Density of the liquid phase, i.e. the one that the vof variable refers to. Recommended value for water is near 1000.

.. input_param:: MultiPhase.density_fluid2

   **type:** Real, optional, default = 1.

   Density of the gas phase, i.e. the one that the vof variable does not refer to. Recommended value for air is near 1.

.. input_param:: MultiPhase.verbose

   **type:** Integer, optional, default = 0

   Verbosity of the MultiPhase physics output. Verbosity greater than 0 will output conservation statistics, recording the
   total water volume and air volume in the domain, the differences between these and the initial volumes, and the difference
   between the initial momentum and the current momentum. These quantities can be used to confirm conservation properties
   in periodic cases without source terms.
