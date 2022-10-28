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
