Theory Manual
================================

Governing equations
------------------------------------

Conservation of fluid mass:

.. math:: \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho U)  = 0

Conservation of fluid momentum:

.. math:: \frac{ \partial (\rho U)}{\partial t} 
   + \nabla \cdot (\rho U U) + \nabla p = \nabla \cdot \tau + \rho g

Incompressibility constraint:

.. math:: \nabla \cdot U = 0

Tracer(s) advection:

.. math:: \frac{\partial \rho s}{\partial t} + \nabla \cdot (\rho U s)  = 0

Discretization 
-------------------------------------
The numerical methdology used to solve the partial differential equations (PDEs)
within AMR-Wind is documented in `Almgren et al. (JCP 1998)
<https://ccse.lbl.gov/Publications/almgren/abchw.pdf>`_.

Time Step -- MOL
~~~~~~~~~~~~~~~~

In the predictor

-  Define :math:`U^{MAC,n}`, the face-centered (staggered) MAC velocity which is used for advection, using :math:`U^n`

-  Define an approximation to the new-time state, :math:`(\rho U)^{\ast}` by setting 

.. math:: (\rho U)^{\ast} &= (\rho U)^n -  
           \Delta t \left( \nabla \cdot (\rho U^{MAC} U) + \nabla {p}^{n-1/2} \right) \\ &+ 
           \Delta t \left( \nabla \cdot \tau^n + \sum_p \beta_p (V_p - {U}^{\ast}) + \rho g \right)

-  Project :math:`U^{\ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t} 
          U^{\ast}+ \frac{1}{\rho} \nabla {p}^{n-1/2} \right)

then defining

.. math:: U^{\ast \ast} = U^{\ast} - \frac{\Delta t}{\rho} \nabla \phi

and 

.. math:: {p}^{n+1/2, \ast} = \phi


In the corrector

-  Define :math:`U^{MAC,\ast \ast}` at the "new" time using :math:`U^{\ast \ast}`

-  Define a new approximation to the new-time state, :math:`(\rho U)^{\ast \ast \ast}` by setting  

.. math:: (\rho U)^{\ast \ast \ast} &= (\rho U)^n - \frac{\Delta t}{2} \left( \nabla \cdot (\rho U^{MAC} U)^n + \nabla \cdot (\rho U^{MAC} U)^{\ast \ast}\right) + \\ &+ \frac{\Delta t}{2} \left( \nabla \cdot \tau^n + \nabla \cdot \tau^{\ast \ast \ast} \right) + \Delta t \left( - \nabla {p}^{n+1/2,\ast} + \sum_p \beta_p (V_p - {U}^{\ast \ast \ast}) + \rho g \right)

-  Project :math:`U^{\ast \ast \ast}` by solving

.. math:: \nabla \cdot \frac{1}{\rho} \nabla \phi = \nabla \cdot \left( \frac{1}{\Delta t} U^{\ast \ast \ast} + \frac{1}{\rho} \nabla {p}^{n+1/2,\ast} \right)

then defining

.. math:: U^{n+1} = U^{\ast \ast \ast} - \frac{\Delta t}{\rho} \nabla \phi

and 

.. math:: {p}^{n+1/2} = \phi

Time Step -- Godunov
~~~~~~~~~~~~~~~~~~~~

-  Define :math:`U^{MAC,n+1/2}`, the MAC velocity which is used for advection. This velocity is interpolated to the cell faces and extrapolated forward in time using a Taylor expansion. Then it is projected to form a divergence-free velocity field.

.. math:: u_f^{n+1/2} = u^n \pm \frac{\Delta x}{2}\frac{\partial u}{\partial x} &-
          \frac{\Delta t}{2}\left(u^n \frac{\partial u}{\partial x}
          + v^n \frac{\partial u}{\partial y} + w^n \frac{\partial u}{\partial z}\right)
          \\ &+
          \frac{\Delta t}{2}\left(g + \frac{1}{\rho^n}\left(
          -\frac{\partial p^{n-1/2}}{\partial x} + \mu\nabla^2u^n\right) \right)

.. math:: \nabla \cdot \frac{\nabla \phi}{\rho^n} = \nabla \cdot \boldsymbol{u}_f^{n+1/2}

.. math:: \boldsymbol{U}^{MAC,n+1/2} = \boldsymbol{u}_f^{n+1/2} - \frac{\nabla \phi}{\rho^n}

-  Time discretization of momentum governing equation

  - The time step goes from :math:`n` to :math:`n+1` with the right-hand-side at :math:`n+1/2`.
  - The time discretization of :math:`\boldsymbol{\tau}` depends on the method chosen to compute diffusion.

.. math:: (\rho \boldsymbol{U})^{n+1} = (\rho \boldsymbol{U})^n &-
           \Delta t \left( \nabla \cdot (\rho \boldsymbol{U} \otimes \boldsymbol{U}^{MAC})^{n+1/2}
           + \nabla {p}^{n+1/2} \right) \\ &+
           \Delta t \left( \nabla \cdot \boldsymbol{\tau} + \rho^{n+1/2} g \right)

-  Partition the discretized equation into two steps, the Predictor and Applying the Projection. The predicted state, :math:`\ast`, is an approximation to the new-time state, :math:`n+1`.

.. math:: (\rho \boldsymbol{U})^{\ast} = (\rho \boldsymbol{U})^n &-
           \Delta t \left( \nabla \cdot (\rho \boldsymbol{U} \otimes \boldsymbol{U}^{MAC})^{n+1/2}
           + \nabla {p}^{n-1/2} \right) \\ &+
           \Delta t \left( \nabla \cdot \boldsymbol{\tau} + \rho^{n+1/2} g \right)

.. math:: (\rho\boldsymbol{U})^{n+1} = (\rho\boldsymbol{U})^{\ast} + \Delta t \nabla p^{n-1/2} - \Delta t \nabla p^{n+1/2}

- In order to calculate :math:`\rho \boldsymbol{U}^{n+1/2}` within the advection term, the momentum must be interpolated to the faces and extrapolated forward in time a half step, similar to the face velocities involved in the MAC projection. To do this for the momentum, the routine uses the momentum, pressure gradients, source terms, and diffusion terms at :math:`n`, as well as :math:`\boldsymbol{U}^{MAC,n+1/2}`.

- In the case of variable density single-phase or multiphase simulations, the density at :math:`n+1` is found using separate scalar equations, which are solved during the predictor step. Because density has no projection step,

.. math:: \rho^{n+1} = \rho^{\ast}.

Therefore, the equation that applies the new pressure gradient becomes

.. math:: \boldsymbol{U}^{n+1} = \boldsymbol{U}^{\ast} + \frac{1}{\rho^{n+1}}\left(\Delta t \nabla p^{n-1/2} - \Delta t \nabla p^{n+1/2}\right)

- The pressure gradient at :math:`n+1/2` is found by solving the projection

.. math:: \nabla \cdot \frac{1}{\rho^{n+1}} \nabla p^{n+1/2} = \nabla \cdot \left( \frac{1}{\Delta t}
          \boldsymbol{U}^{\ast}+ \frac{1}{\rho^{n+1}} \nabla {p}^{n-1/2} \right)

Solving physics on a stretched AMReX mesh
------------------------------------------

Often for simulations involving walls, (e.g., channel flows, complex terrains etc.) it is desirable to have finer mesh near the wall which gradually coarsens only in the wall-normal direction. Consequently, modifications within AMR-Wind are underway to support solving a non-uniformly spaced mesh while still using most of AMReX's machinery directed at uniformly-spaced Cartesian meshes. The governing equations solved on a non-uniform stretched mesh are further explained below -

.. toctree::
   :glob:
   :maxdepth: 2
   
   mapping.rst

Multiphase flow modelling
------------------------------------

The level-set (LS) method was introduced by \cite{OsherSethian1988} to simulate the motion of a surface with curvature dependent speed. In the LS method, the two-phase interface is represented implicitly by a signed distance function (also called level-set function)


.. math::

   \phi ( \mathbf{x},t ) = \begin{cases} 
      d , & \text{in water} \\
      0,  & \text{on surface} \\
      -d, & \text{in air.}
   \end{cases}

where :math:`d` is the distance from point :math:`\mathbf{x}` to the interface.

.. math:: \frac{\partial \phi}{\partial t} + \mathbf{u} \nabla \phi = 0 

where :math:`\mathbf{u}=(u_x,u_y,u_z)` is the velocity vector.
Re-initialization of the level-set

.. math:: \frac{\partial \phi}{\partial \tau} = \text{sgn}(\phi^0)(1-|\nabla \phi|)

where :math:`\phi^0=\phi(x,0)` represents the location of the interface. 

Source terms
------------------------------------

Gravity Forcing
~~~~~~~~~~~~~~~~

The implementation of this source term allows the user to choose the full gravity term (:math:`\rho g`) or a perturbational form (:math:`(\rho - \rho_0) g`). By default, the full term is used, but the perturbational form can be turned on by adding ``ICNS.use_perturb_pressure = true`` to the input file.

The reference density (:math:`\rho_0`) is defined as ``1.0`` by default, can be defined as a constant through the input argument, ``incflo::density``, or can be defined as a spatially varying field within the flow setup (see physics/multiphase/Multiphase.cpp).

Using the perturbational form implies that the hydrostatic pressure is removed from the pressure variable, including its output. This means that the solution to the Poisson equation is actually the perturbational pressure, :math:`p'`, not :math:`p`. If the full pressure, :math:`p`, is desired for analysis or postprocessing purposes, the hydrostatic pressure can be added back to the pressure field via the input argument ``ICNS.reconstruct_true_pressure = true``. In order for this to operate in the code, the reference pressure field must be defined for the specific flow case being run. 

- An example of this is in physics/multiphase/Multiphase.cpp. To construct the reference pressure field, the reference gravity term must be integrated. This particular example assumes that the reference density only varies in z (or is constant), gravity acts only in z, and the hydrostatic pressure at zhi is equal to 0. 

- In mathematical form, the derivation and calculation of the full pressure is as follows:

.. math:: \nabla p = \nabla p' + \rho_0 \boldsymbol{g}

- assume :math:`\boldsymbol{g} = g\hat{k}` and :math:`\frac{dp_0}{dz} = g\hat{k}`

.. math:: p = p' + \int_{z_{min}}^z \rho_0 g dz + p(z = z_{min}) 

- reframe in reference to the top boundary, and assume :math:`p(z = z_{max}) = 0`
   
.. math:: p = p' - \int_z^{z_{max}} \rho_0 g dz + p(z = z_{max}) = p' - \int_z^{z_{max}} \rho_0 g dz


Turbulence Models
-----------------------------

LES models for subgrid scales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Smagorinsky model
^^^^^^^^^^^^^^^^^

Simple eddy viscosity model, the dissipation is calculated using the
resolved strain rate tensor and the grid resolution as

.. math::

   \begin{aligned}
       \tau_{ij} &= -2 \nu_t \widetilde{S}_{ij} \\
       \nu_t &= C_s^2 \Delta^2 (2 \langle S_{ij} S_{ij} \rangle)^{\frac{1}{2}}
   \end{aligned}


AMDNoTherm model
^^^^^^^^^^^^^^^^^
This is the implementation of the base AMD model, useful for flows without a temperature field.

The eddy viscosity is calculated using an anisotropic derivative with a
different filter width in each direction

.. math::

   \begin{aligned}
       \hat{\partial}_i &= \sqrt{C} \delta_i \partial_i \textrm{ for } i=1,2,3 \\
       C &= 1/3, \textrm{ Poincare coefficient for } 2^{nd} \textrm{ order gradient} \\
       \delta_i &= \textrm{Filter width along dimension } i \textrm{ for anisotropic grids}
   \end{aligned}

The anisotropic derivative is used to define the eddy viscosity as

.. math::

   \begin{aligned}
       \tau_{ij} &= -2 \nu_t \widetilde{S}_{ij} \\
       \nu_t &= \frac{- (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{u}_j) \widetilde{S}_{ij}}{ (\partial_l \widetilde{u}_m) (\partial_l \widetilde{u}_m) }
   \end{aligned}


AMD model (for ABL)
^^^^^^^^^^^^^^^^^^^

The eddy viscosity is calculated using an anisotropic derivative with a
different filter width in each direction

.. math::

   \begin{aligned}
       \hat{\partial}_i &= \sqrt{C} \delta_i \partial_i \textrm{ for } i=1,2,3 \\
       C &= 1/3 \textrm{ Poincare coefficient for } 2^{nd} \textrm{ order gradient} \\
       \delta_i &= \textrm{Filter width along dimension } i \textrm{ for anisotropic grids}\\
       \beta &= g/\Theta_0 \textrm{ Gravity constant over reference temperature}
   \end{aligned}

The anisotropic derivative is used to define the eddy viscosity as

.. math::

   \begin{aligned}
       \tau_{ij} &= -2 \nu_t \widetilde{S}_{ij} \\
       \nu_t &= \frac{- (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{u}_j) \widetilde{S}_{ij} +  \beta (\hat{\partial}_k \widetilde{w}) (\hat{\partial}_k (\widetilde{\Theta} - \langle {\widetilde{\Theta}} \rangle) )  }{ (\partial_l \widetilde{u}_m) (\partial_l \widetilde{u}_m) } \\
       \tau_{\theta j} &= -2 D_e \frac{\partial \widetilde{\Theta}}{\partial x_j} \\
       D_e &= \frac{- (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{\Theta}) \partial_i \widetilde{\Theta} }{(\partial_l \widetilde{\Theta}) (\partial_l \widetilde{\Theta})}
   \end{aligned}

- **Implementation details:**


For ease of implementation, each part of :math:`\nu_t` and :math:`D_e`
is expanded in this subsection, these are used in ``AMD.h`` within
functions ``amd_muvel`` and ``amd_thermal_diff``.

**Terms for** :math:`\nu_t` **in** ``amd_muvel``


#. :math:`(\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{u}_j) \widetilde{S}_{ij}`
   is a contraction of 2 symmetric tensors, :math:`(\hat{\partial}_k
   \widetilde{u}_i) (\hat{\partial}_k \widetilde{u}_j)` and
   :math:`\widetilde{S}_{ij}`, therefore we get 6 unique terms, 3
   diagonals and 3 off-diagonals. The diagonal terms get a factor of 2
   from :math:`\widetilde{S}_{ij}` and the off-diagonal terms get a
   factor of 2 from symmetry. This term is ``num_shear``.

   .. math::

      \begin{aligned}
          \begin{split}
          (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{u}_j) \widetilde{S}_{ij} = \\
          2*C* [\\
          {}&u_x*(u_x^2 dx^2 + u_y^2 dy^2 + u_z^2 dz^2) +\\
          {}&v_y*(v_x^2 dx^2 + v_y^2 dy^2 + v_z^2 dz^2) +\\
          {}&w_z*(w_x^2 dx^2 + w_y^2 dy^2 + w_z^2 dz^2) +\\
          {}&(u_y+v_x) * (
          u_x v_x dx^2 +
          u_y v_y dy^2 +
          u_z v_z dz^2
          ) +\\
          {}&(u_z+w_x) * (
          u_x w_x dx^2 +
          u_y w_y dy^2 +
          u_z w_z dz^2
          ) +\\
          {}&(v_z+w_y) * (
          v_x w_x dx^2 +
          v_y w_y dy^2 +
          v_z w_z dz^2
          )\\
          ]
          \end{split}
      \end{aligned}

#. :math:`\beta (\hat{\partial}_k \widetilde{w}) (\hat{\partial}_k (\widetilde{\Theta} - \langle {\widetilde{\Theta}} \rangle) )`
   is implemented as ``num_buoy``

   .. math::

      \beta (\hat{\partial}_k \widetilde{w}) (\hat{\partial}_k (\widetilde{\Theta} - \langle {\widetilde{\Theta}} \rangle) ) =
      \beta* C* ( w_x \Theta_x dx^2 + w_y \Theta_y dy^2 + w_z \Theta_z dz^2)

#. :math:`(\partial_l \widetilde{u}_m) (\partial_l \widetilde{u}_m)` is
   double contraction of rank 2 tensors, and has 9 unique terms. This is
   implemented as ``denom``

   .. math::

      \begin{aligned}
          \begin{split}
              (\partial_l \widetilde{u}_m) (\partial_l \widetilde{u}_m) = \\
              {}& u_x u_x + u_y u_y + u_z u_z \\
              {}& v_x v_x + v_y v_y + v_z v_z \\
              {}& w_x w_x + w_y w_y + w_z w_z
          \end{split}
      \end{aligned}

**Terms for** :math:`D_e` **in** ``amd_thermal_diff``

#. :math:`(\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{\Theta}) \partial_i \widetilde{\Theta}`
   is a double contraction, has 9 unique terms. This is implemented as
   ``num``

   .. math::

      \begin{aligned}
          \begin{split}
          (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{\Theta}) \partial_i \widetilde{\Theta} = \\
          C*[ \\
          {}& (u_x \Theta_x dx^2 + u_y \Theta_y dy^2 + u_z \Theta_z dz^2)*\Theta_x +\\
          {}& (v_x \Theta_x dx^2 + v_y \Theta_y dy^2 + v_z \Theta_z dz^2)*\Theta_y +\\
          {}& (w_x \Theta_x dx^2 + w_y \Theta_y dy^2 + w_z \Theta_z dz^2)*\Theta_z \\
          ]
          \end{split}
      \end{aligned}

#. :math:`(\partial_l \widetilde{\Theta}) (\partial_l \widetilde{\Theta}) = \Theta_x^2 + \Theta_y^2 + \Theta_z^2`
   is implemented as ``denom``

- **Unit tests**

There is a simple unit test for both :math:`\nu_t` and :math:`D_e` in
``unit_tests/turbulence/test_turbulence_LES.cpp`` under
``test_AMD_setup_calc``.

Wall models
-----------
The wall models descibed in this section are implemented in ``AMR-wind`` for
running wall-bounded flows (non-ABL cases).

Log-law wall model
~~~~~~~~~~~~~~~~~~

This wall model computes the local :math:`u_\tau` from the velocity at
the first grid cell, and uses this to compute the shear stress, which is
then used as a boundary condition.

The log law:

.. math:: u_{\mathrm{mag}} = u_\tau \left(\frac{1}{\kappa}\log\left(\frac{u_\tau z}{\nu}\right) + B\right). \label{eq:loglaw}

Given a horizontal velocity magnitude
:math:`u_{\mathrm{mag}} = \sqrt{u^2 + v^2}` at
:math:`z = z_{\mathrm{ref}}`, :math:`u_\tau` can be computed using a
non-linear solve to satisfy `[eq:loglaw] <#eq:loglaw>`__.

In ``AMR-wind`` Newton-Raphson iterations are used with a convergence
criterion of :math:`\lvert u_\tau^{n+1} - u_\tau^n \rvert < 10^{-5}`.
For this, derivative of
:math:`\frac{\partial u_{\mathrm{mag}}}{\partial {u_\tau}}` is used,

.. math:: \frac{\partial u_{\mathrm{mag}}}{\partial {u_\tau}} = \left(\frac{1}{\kappa}\left(1+\log\left(\frac{u_\tau z_{\mathrm{ref}}}{\nu}\right)\right) + B\right)

.. math:: u_\tau^{n+1} = u_\tau^{n} - \left(u_\tau^n \left(\frac{1}{\kappa}\log\left(\frac{u_\tau^n z_{\mathrm{ref}}}{\nu}\right) + B\right) - u_{\mathrm{mag}}\right)/\frac{\partial u_{\mathrm{mag}}}{\partial {u_\tau}}.

Finally, the shear stress is calculated as,

.. math::

   \begin{aligned}
       \tau_{xz} &= u_\tau^2 \frac{u}{u_\mathrm{mag}} \\
       \tau_{yz} &= u_\tau^2 \frac{v}{u_\mathrm{mag}}
   \end{aligned}

Constant stress model
~~~~~~~~~~~~~~~~~~~~~

NOTE: This wall model will be ill-posed unless combined with a Dirichlet
boundary condition on the other wall, :math:`\langle u \rangle` can
drift by a constant otherwise.

This is a trivial wall model, where the shear stresses are specified as
constants. For a pressure gradient driven channel,

.. math::

   \begin{aligned}
       u_\tau^2 &= -\frac{\mathrm{d} P}{\mathrm{d} x} \\
       \tau_{xz} &= u_\tau^2 \\
       \tau_{yz} &= 0
   \end{aligned}

Schumann model
~~~~~~~~~~~~~~

NOTE: This wall model will be ill-posed unless combined with a Dirichlet
boundary condition on the other wall, :math:`\langle u \rangle` can
drift by a constant otherwise.

This model is a modified version of the constant stress model, where the
fluctuations from a reference height :math:`z_\mathrm{ref}` are used to
add fluctuations in the shear stress.

.. math::

   \begin{aligned}
       u_\tau^2 &= -\frac{\mathrm{d} P}{\mathrm{d} x} \\
       \tau_{xz} &= u_\tau^2 \frac{u}{\langle u_\mathrm{mag} \rangle} \\
       \tau_{yz} &= u_\tau^2 \frac{v}{\langle u_\mathrm{mag} \rangle}
   \end{aligned}

where, :math:`\langle u_\mathrm{mag} \rangle` is the planar average of
:math:`u_{\mathrm{mag}} = \sqrt{u^2 + v^2}` at :math:`z_\mathrm{ref}`.

Symmetric wall boundary
~~~~~~~~~~~~~~~~~~~~~~~

This is a boundary condition to for flows with a symmetry across the
z direction (example: *half-channel* simulations) at the centerline.

.. math::

   \begin{aligned}
       \tau_{xz} &= 0 \\
       \tau_{yz} &= 0 \\
       w &= 0
   \end{aligned}

Navigating source code
------------------------

``AMR-Wind`` is built on top of `AMReX library
<https://amrex-codes.github.io/amrex/docs_html/>`_. Users are strongly
recommended to read through the AMReX documentation and understand the basic
AMReX concepts before jumping into the AMR-Wind source code.

The `Basics section
<https://amrex-codes.github.io/amrex/docs_html/Basics_Chapter.html>`_ provides a
thorough overview of the basic datastructures and ways to interact with these
structures. The `GPU section
<https://amrex-codes.github.io/amrex/docs_html/GPU_Chapter.html>`_ provides an
overview of the AMReX GPU strategy and the higher-level functions (e.g.,
``parallel-for`` abstractions) available to write GPU-ready code within
AMR-Wind. The `Linear Solvers section
<https://amrex-codes.github.io/amrex/docs_html/LinearSolvers_Chapter.html>`_
provides an overview of the multi-level multigrid (MLMG) solvers used to solve
the various linear systems within AMR-Wind.
