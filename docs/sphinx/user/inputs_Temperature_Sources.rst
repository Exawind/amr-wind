.. _inputs_temperature_sources:
   
Section: Temperature Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: temperature.source_terms

   **type:** String(s), optional
   
   Activates source terms for the energy equations. These strings can be 
   entered in any order with a space between
   each. Please consult the :doc:`../doxygen/html/index` for a
   comprehensive list of all energy source terms available. Note that the
   following input arguments specific to each source term will only be active
   if the corresponding source term (the root name) is listed in 
   :input_param:`temperature.source_terms`.

.. input_param:: DragTempForcing.drag_coefficient

   **type:** Real, optional

   This value specifies the coefficient for the forcing term in the immersed boundary forcing method. It is currently
   recommended to use the default value to avoid initial numerical stability. 

.. input_param:: DragTempForcing.bc_forcing_time_factor

   **type:** Real, optional, default = 5.0

   This value modifies the time scale of the BC forcing component of DragTempForcing relative to
   the time step size.


The following arguments are influential when ``ImmersedDragTempForcing`` is included in
``Temperature.source_terms``. This is the temperature counterpart of
``ImmersedDragForcing`` and requires the :ref:`ImmersedTerrain <inputs_immersedterrain>`
physics. Inside the terrain the temperature relaxes toward the soil temperature at the rate
:math:`\beta C_d / \Delta z`, integrated exactly over the time step. When
:input_param:`turbulence.model` is not ``Laminar`` a Monin-Obukhov heat-flux wall model is
applied in surface cells with weight :math:`1 - \beta`, on the same wall patches and with the
same friction velocity as the momentum source: the wall-model geometry, relaxation time scale
and minimum roughness are read from :input_param:`ImmersedDragForcing.wall_model`,
:input_param:`ImmersedDragForcing.bc_forcing_time_scale`,
:input_param:`ImmersedDragForcing.bc_forcing_time_factor` and
:input_param:`ImmersedDragForcing.minimum_z0`. The forcing is the surface flux divergence
:math:`-u_* \theta_* / \Delta_n` plus a relaxation toward the log-law temperature.

.. input_param:: ImmersedDragTempForcing.drag_coefficient

   **type:** Real, optional, default = 10.0

   Coefficient :math:`C_d` of the interior relaxation; the rate in a fully solid cell is
   :math:`C_d / \Delta z`. Unlike ``DragTempForcing`` the rate does not depend on the local
   velocity, following the constant velocity scale of Muñoz-Esparza and coworkers (2020).

.. input_param:: ImmersedDragTempForcing.soil_temperature

   **type:** Real, optional, default = 300.0

   Temperature the terrain interior relaxes toward. Also the surface temperature when
   ``surface_condition = surface_temperature``.

.. input_param:: ImmersedDragTempForcing.surface_condition

   **type:** String, optional, default = ``obukhov_length``

   How the temperature scale :math:`\theta_*` of the wall model is obtained.

   - ``obukhov_length``: :math:`\theta_* = \theta u_*^2 / (\kappa g L)` from the single
     :input_param:`ABL.monin_obukhov_length`; the surface temperature is inferred from the
     reference cell. Reproduces ``DragTempForcing``.
   - ``surface_temperature``: the surface is held at ``soil_temperature`` and
     :math:`\theta_* = \kappa (\theta_\mathrm{ref} - \theta_s) / \phi_h(d_2)`, so the interior
     and the surface use the same temperature.
   - ``heat_flux``: :math:`\theta_* = -q_s / u_*` from
     :input_param:`ImmersedDragTempForcing.surface_heat_flux`.

.. input_param:: ImmersedDragTempForcing.surface_heat_flux

   **type:** Real, optional, default = 0.0

   Kinematic surface heat flux :math:`q_s = \overline{w'\theta'}_s` (K m/s, positive heats the
   air) for ``surface_condition = heat_flux``.

.. input_param:: ImmersedDragTempForcing.force_laminar

   **type:** Boolean, optional, default = false

   Skip the wall model even when a turbulence model is active.

The following list of inputs are used with the `Temperature.source_terms = PerturbationForcing` option to add perturbation to the 
temperature field to generate flow structures for LES when the inflow data is coarse or uniform flow condition. Not 
recommended for use with RANS models. 

.. input_param:: PerturbationForcing.start

   **type:** Real, mandatory

   Start location of the perturbation box 

.. input_param:: PerturbationForcing.end

   **type:** Real, mandatory

   End location of the perturbation box 

..  input_param:: PerturbationForcing.pert_amplitude

   **type:** Real, optional 

   Amplitude of temperature perturbation 

..  input_param:: PerturbationForcing.time_steps 

   **type:** Real, optional 

   Separation time between applying perturbations. A high value may dampen the flow structures 
   and a small value may cause numerical instability. 