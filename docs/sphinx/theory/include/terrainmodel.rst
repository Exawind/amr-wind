.. _terrainmodel:

Terrain Model
--------------

An immersed boundary forcing method (IBFM) is used to represent the terrain. In this method,
the effect of the terrain is modeled using a forcing term in the momentum and energy equation.
Two implementations are available: the original ``TerrainDrag`` physics with a binary blanking
of cells, and the ``ImmersedTerrain`` physics with a partial terrain fraction per cell. Both
follow the immersed body force method of
`Muñoz‐Esparza, Domingo, et al. (JAMS 2020) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002141>`_,
which prescribes only the forcing inside the body; the wall functions, stability corrections and
time integration described below are additions made in this code.

Binary blanking (``TerrainDrag``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The forcing term in the momentum equation is given by:

.. math::

   F_i = - \beta C_d u_i | u_i |

Here :math:`\beta` is the volume fraction of the cell covered by terrain, :math:`C_d` is a drag
term  and :math:`u_i` is the wind speed. In ``TerrainDrag`` the volume fraction is
computed as a 0 or 1 using a simple nearest cell algorithm at each grid level, which turns slopes
into a staircase. The calculation of the drag coefficient term and the forcing term for the energy
equation can be found in the reference above.

The original formulation is designed for low Reynolds number cases and does not include a
method for applying a wall function. We propose the use of a forcing function to include
the wall effects.

First, compute the friction velocity from the cell above the terrain-adjacent cell, whose
center is :math:`1.5\,\Delta z` above the wall. The velocity of the terrain-adjacent cell itself
is not used because that is the cell being forced:

.. math::

   u_*= |u_h[k+1]| \frac {\kappa}{\log [1.5 \Delta z/z_0] - \psi_m(1.5 \Delta z / L)}

The expected wind speed at cell k, whose center is :math:`0.5\,\Delta z` above the wall, is

.. math::

   |u_n|= \frac{u_*}{\kappa} \left[ \log (0.5 \Delta z/z_0) - \psi_m(0.5 \Delta z / L) \right]

with :math:`\psi_m` the Monin-Obukhov stability function for a single prescribed Obukhov
length :math:`L` (neutral when it is not specified). The forcing term is computed as

.. math::

   F_i= - \frac {|u[k]| \hat{c} - |u_n|\hat{l}} {\tau}

Here :math:`\hat{c}=(1,1,1)` is the existing normal vector from the grid and :math:`\hat{l}=(ux,uy,0)/|u_n|` is the value
from the log law.

.. image:: ./images/terrain_normal.png
   :align: center
   :width: 30%

Partial terrain fraction (``ImmersedTerrain``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ImmersedTerrain`` replaces the binary blanking with the fraction :math:`\beta \in [0,1]` of each
cell occupied by terrain (the fraction of the cell column below the terrain height, or optionally
a smooth function of the signed distance to the surface). It also stores an integer mask
(0 fluid, 1 solid, 2 surface), the terrain height and its slopes
:math:`\partial h/\partial x, \partial h/\partial y`, and the roughness length. Surface cells are the
partially filled cells and the fluid cells sharing a face with a mostly solid cell on any of the six
sides, so that the vertical faces of steep terrain and buildings are treated as walls as well.
The mask can be used directly by ``FieldRefinement`` to refine the surface band.

**Immersed drag.** The drag inside the terrain is applied by ``ImmersedDragForcing`` as a linear
relaxation toward zero velocity at the rate :math:`C = \beta C_d/\Delta z`. Because the term is
linear in :math:`u`, it is integrated exactly over the time step,

.. math::

   \frac{\partial u_i}{\partial t} = - C_\mathrm{eff}\, u_i, \qquad
   C_\mathrm{eff} = \frac{1 - e^{-C \Delta t}}{\Delta t} \le \min\left(C, \frac{1}{\Delta t}\right),

which is stable and monotone for any :math:`C \Delta t` and removes the need for drag limiters.
The explicit form :math:`-C u_i` overshoots at cold start whenever :math:`C \Delta t > 2`, which
happens routinely on the first CFL-limited step.

**Implicit projection.** With the drag applied explicitly, the velocity inside the body is
reset to zero by the source term and then re-created by the pressure projection, which treats the
body as fluid and adds :math:`\Delta t \nabla p / \rho` back every step. Under a CFL time step this
residual is proportional to :math:`\Delta x`, so the approach to zero inside the body is first order
regardless of how the interface is represented. The option ``ImmersedTerrain.implicit_projection``
applies the drag implicitly together with the pressure,

.. math::

   u^{n+1} = \frac{u^{**} - \Delta t \nabla p / \rho}{1 + \beta C \Delta t},
   \qquad \nabla \cdot u^{n+1} = 0
   \;\Rightarrow\;
   \nabla \cdot \left( \frac{\Delta t}{\rho (1 + \beta C \Delta t)} \nabla p \right)
   = \nabla \cdot \frac{u^{**}}{1 + \beta C \Delta t},

that is, the terrain behaves as a fluid of density :math:`\rho (1 + \beta C \Delta t)` in the
nodal and MAC projections. The residual inside the body becomes independent of the time step and
equal to :math:`\nabla p / (\rho C)`; increasing :math:`C_d` in proportion to :math:`1/\Delta z`
then gives second-order convergence. On the laminar immersed box case the fitted orders of the
mean speed inside the body were 1.3 (binary), 1.1 (partial fraction, explicit drag), 1.0
(implicit projection, fixed :math:`C_d`) and 1.9 (implicit projection, :math:`C_d \propto 1/\Delta z`).

**Wall model.** In surface cells, weighted by :math:`1 - \beta` and only when a turbulence model
is active, a Monin-Obukhov wall model is applied on one or more wall patches. Each patch has a wall
distance :math:`d_1` of the cell center, a distance :math:`d_2` of a reference cell whose velocity
is trusted, a cell width :math:`\Delta_n` along the wall normal and a unit normal :math:`\hat n`
pointing into the fluid. Three ways of building the patches are available through
``ImmersedDragForcing.wall_model``:

* ``cell_offset``: one patch per face touching a mostly solid neighbor, reference cell opposite
  the wall, :math:`d_1 = 0.5 \Delta_f`, :math:`d_2 = 1.5 \Delta_f`;
* ``terrain_height``: as above, but the bottom face uses the true height above the terrain,
  :math:`d_1 = z_k - h`, :math:`d_2 = z_{k+1} - h`;
* ``surface_normal``: a single patch along
  :math:`\hat n = (-h_x, -h_y, 1)/\sqrt{1 + h_x^2 + h_y^2}` with :math:`d_1 = (z_k - h)\, n_z`, the
  reference cell one normal cell-width away, and the velocity split into wall-normal and tangential
  parts so that the stress acts on :math:`u, v, w` together.

With :math:`u_t` the tangential velocity and :math:`\phi_m(d) = \log(d/z_0) - \psi_m(d/L)`,

.. math::

   u_* = \frac{\kappa |u_{t,\mathrm{ref}}|}{\phi_m(d_2)}, \qquad
   |u_{t,1}| = \frac{u_*}{\kappa} \phi_m(d_1), \qquad
   F_t = - \frac{u_*^2}{\Delta_n} \hat e_t - C_{bc} \left( u_t - |u_{t,1}| \hat e_t \right),

where :math:`\hat e_t` is the direction of the tangential reference velocity. Only the
tangential velocity is relaxed; the wall-normal component is left to the projection. The
relaxation rate :math:`C_{bc}` uses the same exact integration as the drag, with time scale
:math:`\max(\tau_f \Delta t,\, d_1/u_*)` by default, so that the wall model becomes independent of
the time step once the flow-based scale exceeds the floor. Contributions from several patches are
averaged with the neighbor fractions as weights.

**Energy equation.** ``ImmersedDragTempForcing`` relaxes the temperature inside the terrain toward
the soil temperature at the rate :math:`\beta C_d/\Delta z`, without dependence on the local
velocity, and applies a heat-flux wall model on the same patches and with the same friction
velocity as the momentum source:

.. math::

   \theta_\mathrm{target} = \theta_s + \frac{\theta_*}{\kappa} \phi_h(d_1), \qquad
   F_\theta = \frac{-u_* \theta_*}{\Delta_n} - C_{bc} \left( \theta - \theta_\mathrm{target} \right),

with :math:`\theta_*` obtained from a prescribed Obukhov length
(:math:`\theta_* = \theta u_*^2 / (\kappa g L)`), a prescribed surface temperature
(:math:`\theta_* = \kappa (\theta_\mathrm{ref} - \theta_s)/\phi_h(d_2)`), or a prescribed surface
heat flux (:math:`\theta_* = -q_s/u_*`).

**Implicit drag in the diffusion solve.** With ``implicit_projection`` the terrain cells must
also be held at rest inside the implicit diffusion solve: the operator uses the same effective
density :math:`\rho(1 + \beta C \Delta t)` as the coefficient of the time-derivative term (with the
right-hand side kept at the plain density), so that momentum diffusing into the terrain is damped
there rather than accumulated and removed by the projection. Without this the terrain cells float
to a fraction of the fluid velocity during each solve, the interface flux is reduced, the effective
wall sits about one cell too low, and the coupled scheme becomes unstable when
:math:`\nu \Delta t / \Delta z^2` exceeds about 10.

**Diffusive flux at the interface.** The diffusion operator evaluates the viscous flux at a
fluid/solid face as :math:`\mu_\mathrm{eff}(u_k - 0)/\Delta_f`, since the interior is at rest.
This is a wall stress with the wrong length scale: the flux to a no-slip wall at distance
:math:`d_1` is :math:`\mu u_k/d_1`, and with a turbulence model the stress is already supplied
by the wall model, so the SGS flux across the interface counts it twice. The option
``ImmersedTerrain.interface_diffusion`` modifies the face coefficients: ``block`` multiplies
them by :math:`\min(1-\beta_L, 1-\beta_R)`, removing the flux across the interface so the wall
model alone acts (turbulent pathway); ``no_slip`` multiplies fluid/solid faces by
:math:`\Delta_f/d_1` with :math:`d_1 = z_k - h` on the bottom face, which reproduces the
no-slip flux at the true wall position (laminar pathway with finite viscosity). With the
default ``none`` the molecular viscosity should be kept negligible in laminar test cases.

**SGS viscosity at the interface (Kosovic model).** The ``Kosovic`` model carries its own
interface treatment: the SGS viscosity is zero inside the terrain and, in the first fluid
cells, replaced by :math:`2 \rho u_*^2 / |\partial U_t/\partial n|`, so that the SGS flux across
the interface (face viscosity equal to the mean of the cell and the blanked neighbor) equals the
wall stress :math:`\rho u_*^2` instead of the unrelated value the LES closure would give there.
With ``Kosovic.terrain_model = ImmersedTerrain`` this is evaluated from the terrain fraction and
the surface cells: the viscosity and the non-linear term are weighted by :math:`1 - w_\mathrm{solid}`
and the log-law value is averaged over the wall patches of the surface cell, with the same
reference cell, distance and normal as the wall model (six faces or surface normal). The
default ``TerrainDrag`` keeps the binary blanking and the drag-cell treatment above the terrain.

**Laminar channel verification.** For plane Poiseuille flow with an immersed flat bottom wall
between cell centers and a no-slip top wall, the combination ``implicit_projection``,
``interface_diffusion = no_slip`` and ``drag_weight = center`` converges to the exact parabolic
profile at second order (L2 error falling by a factor of four per refinement over three
successive refinements, independent of the sub-cell wall position), provided the drag coefficient is large
enough that the residual velocity of the terrain cells, of order :math:`1/(1 + C \Delta t)`, stays
below the discretization error, or is increased with resolution. The default ``none`` places the
effective wall at the solid cell center (first order) and ``drag_weight = fraction`` damps partial
cells whose center lies in the fluid (a mesh-alignment-dependent O(1) error in that cell).
