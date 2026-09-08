# Immersed Body Force Method: reference notes for TerrainDrag / DragForcing / DragTempForcing

Source paper (the basis for the terrain and building forcing in this code):

> Muñoz-Esparza, D., Sauer, J. A., Shin, H. H., Sharman, R., Kosović, B., Meech, S.,
> García-Sánchez, C., Steiner, M., Knievel, J., Pinto, J., Swerdlin, S. (2020).
> *Inclusion of Building-Resolving Capabilities Into the FastEddy GPU-LES Model Using an
> Immersed Body Force Method.* J. Adv. Model. Earth Syst., 12, e2020MS002141.
> https://doi.org/10.1029/2020MS002141

The paper itself extends Chan & Leach (2007), J. Appl. Meteor. Climatol. 12, 2127-2146.

Scope: the paper supplies only the blanked-cell momentum and temperature forcing (Section 1).
The code deliberately goes beyond it with engineering fixes: a log-law wall function in the
first fluid cell, Monin-Obukhov stability corrections, drag limiters for start-up numerical
stability, sponge and Rayleigh damping layers, and a moving-terrain path for ocean waves
(Sections 2-4). Differences from the paper listed below are intentional unless flagged.

Code files:

- `src/physics/TerrainDrag.cpp` (builds the blanking / drag / roughness / height / damping fields)
- `src/equation_systems/icns/source_terms/DragForcing.cpp` (momentum source)
- `src/equation_systems/temperature/source_terms/DragTempForcing.cpp` (temperature source)
- Theory page: `docs/sphinx/theory/include/terrainmodel.rst`

---

## 1. What the paper prescribes

### 1.1 Original IBFM (Chan & Leach 2007)

Momentum forcing in every immersed grid point:

    F_ui = - C_d * rho * |u| * u_i          (Eq. 1)

- `C_d` has units of 1/m (it is really `C_d' * A_p`, a drag coefficient times a plant-area
  density, borrowed from canopy models).
- Chan & Leach tested `C_d` = 15, 50, 100 1/m and found 100 1/m reproduced a solid-wall
  CFD reference. Velocity inside the body goes to "essentially zero".

### 1.2 Extension 1: scale-independent drag coefficient

By balancing the body force against the advection term (Eq. 2):

    C_d = alpha_m * beta_r / Delta,   Delta = (dx dy dz)^(1/3),   alpha_m = 1e3

- `beta_r` is the immersed volume fraction of the cell (1 inside, 0 in fluid, 0-1 partial).
- Lower `alpha_m` lets air leak through the body; larger `alpha_m` caused numerical
  instabilities. `alpha_m = 1e3` gives `C_d` = 100-1000 1/m for Delta = 1-10 m, consistent
  with Chan & Leach.
- In practice the paper found no scale correction was needed for Delta = 2-20 m, so the
  operational form is `C_d = alpha_m * max(1/Delta, 1)`; the correction only matters for
  Delta < 1 m (wind-tunnel scale).
- Residual velocities inside bodies: 1e-3 to 1e-4 m/s. Leakage of a passive tracer into
  buildings was < 5 % of the instantaneous total.

### 1.3 Extension 2: thermal (and mass) forcing

To stop leaked air from homogenising the temperature across the body interface:

    F_theta = - C_t * |U_s| * (rho*theta - (rho*theta)_ref)       (Eq. 3)
    F_rho   = - C_t * |U_s| * (rho - rho_ref)                     (Eq. 4, compressible only)

- `U_s` = 1 m/s is a **constant** velocity scale (deliberately not the local velocity, which
  tends to zero inside the body and would remove the forcing).
- `C_t = alpha_t * beta_r / Delta` with `alpha_t = 10`.
- Verified on an isolated 120 m cube at Delta = 5 m: building held at theta_ref = 300 K
  with 1e-3 K mean error; without the term the interior warmed 1.75 K in 2 h.

### 1.4 What the paper explicitly does NOT do

- **No wall model.** The IBFM is a low-Reynolds-number, no-slip approach. Section 4 attributes
  part of the residual wind-speed error to "the fact that the IBFM does not include a wall
  model", and the conclusions state a plan to add a rough-wall model following a canopy
  formulation (Anderson 2013, International Journal for Numerical Methods in Fluids 71, 1588-1608).
- No partial-cell treatment was exercised (beta_r taken as 0/1 in the validations).
- Lateral forcing is by Dirichlet inflow from a precursor plus the cell-perturbation method,
  not by sponge layers.
- Numerics: 5th-order Wicker-Skamarock advection, RK3 in time, Smagorinsky-Lilly SGS.
- Validation: staggered cube array (Castro et al. 2006, Delta = 1 mm) and Joint Urban 2003
  IOP9 downtown Oklahoma City (Delta = 2 m, 20 sonic anemometers, SF6 release).

---

## 2. How the code maps onto the paper

### 2.1 Blanking (`beta_r`) - TerrainDrag.cpp

- `terrain_blank` = 1 when the cell centre is at or below the bilinearly interpolated terrain
  height (or, for OceanWaves without `vof`, when the wave volume fraction >= 0.5). Pure 0/1;
  no partial volume fraction, as in the paper's validations.
- `terrain_drag` = +1 for the first fluid cell above terrain, -1 for a fluid cell directly
  below terrain (overhangs / channel ceilings). The sign gives the offset to the fluid
  neighbour used by the wall function. This marker field has no counterpart in the paper.
- `terrainz0` (roughness), `terrain_height`, `terrain_damping` also have no counterpart.

### 2.2 Momentum forcing inside blanked cells - DragForcing.cpp

    Cd     = drag_coefficient / dz                 (default drag_coefficient = 10)
    CdM    = min( Cd / (|u| + eps), cd_max / scale_factor )
    CdM_m  = CdM * |u|                              (optionally min(., 1/dt))
    src   -= CdM_m * (u_n - u_target_n) * terrain_blank

Differences from the paper:

| Item | Paper | Code |
|---|---|---|
| Length scale | Delta = (dx dy dz)^(1/3) | dz only |
| alpha | 1e3 (momentum) | 10 (`drag_coefficient`) |
| Velocity dependence | quadratic, C_d |u| u | effectively **linear**: CdM_m ~ Cd, so F = -(10/dz) u. The |u| cancels except when the `cd_max` cap is active at very small |u| |
| |u| used | 3-D speed | 3-D speed (same) |
| Target | 0 | 0, or the wave orbital velocity when terrain is waves |
| Small-Delta correction | C_d = alpha max(1/Delta, 1) | `use_original_drag_limiter`: for dz < 1 the /dz scaling is dropped in the cap (and in Cd itself when `is_laminar`) |

The effective forcing is a Rayleigh-type relaxation of u toward the target with time scale
dz / drag_coefficient (e.g. 1 s for dz = 10 m). This is closer in spirit to
Smolarkiewicz et al. (2007), which the paper cites as the equivalent relaxation approach.

### 2.3 Temperature forcing inside blanked cells - DragTempForcing.cpp

    Cd   = min( (drag_coefficient/dz) / (|u| + tiny), 10 / dz )    (default drag_coefficient = 1)
    src -= Cd * (theta - soil_temperature) * terrain_blank

- The cap `10/dz` is exactly the paper's `alpha_t * U_s / Delta` with alpha_t = 10, U_s = 1 m/s
  and Delta -> dz. Inside a body where |u| -> 0 the cap is active, so the paper's constant-
  velocity-scale form is recovered.
- Away from the cap the coefficient is 1/(dz |u|), i.e. it decreases with wind speed. The
  paper argued against any dependence on the local velocity.
- The paper's mass/density term (Eq. 4) has no analogue; the solver is incompressible.
- `theta_ref` is a single `soil_temperature` (default 300 K), not a per-building value.

---

## 3. Code additions beyond the paper: wall function

The paper has no wall model. The code adds one in the `terrain_drag` cells (|marker| = 1),
skipped when `DragForcing.is_laminar = true`.

### 3.1 Momentum (DragForcing.cpp, `viscous_drag_calculations` and the bc forcing block)

With `k_nb = k + terrain_drag(i,j,k)` the fluid neighbour, `dz` the cell height, wall at the
face between the drag cell and the blanked cell:

    u*        = kappa * |u_h[k_nb]| / ( ln(1.5 dz / z0) - psi_m(1.5 dz / L) )
    Dxz, Dyz  = - u*^2 * (u1, v1) / |u_h1| / dz          (stress divergence over the cell)
    |u_n|     = u*/kappa * ( ln(0.5 dz / z0) - psi_m(0.5 dz / L) )
    u_target  = |u_n| * (u2, v2)/|u_h2|                    (direction of the neighbour)
    bc_force  = - (u_target - u1) / (bc_forcing_time_factor * dt)
    src_x,y  -= (Dxz, Dyz) + bc_force
    src_z    -= CdM_m * (w1 - w_target)                    (w damped like a blanked cell)

- `u_h` is horizontal speed; u* uses the neighbour at 1.5 dz, the target is at 0.5 dz.
- `z0` = max(`terrainz0`, `DragForcing.minimum_z0` = 1e-4). Uniform or from a roughness file.
- `bc_forcing_time_factor` (default 5) relaxes the cell toward the log-law value over five
  time steps rather than imposing it directly.
- For waves the velocities are relative to the wave surface velocity, and an optional
  inviscid form drag (`wave_model_inviscid_form_drag`, MOSD of Ayala et al. 2024) is added.

Why 1.5 dz (rationale from H. Gopalan): the wind field in the first cell above the surface
(the drag cell) is not trustworthy, because that is the cell being forced. The cell above it
is well resolved, so its velocity is used to compute u*, and that u* is then used to force
the drag cell toward the log-law value at its own centre (0.5 dz). The neighbour centre sits
1.5 dz above the wall face, hence ln(1.5 dz / z0). The theory page `terrainmodel.rst` writes
this as `ln((z_{k+1}-z_k)/z0)`, which should read 1.5 dz; the page also predates the
stability functions and the bc time factor.

### 3.2 Temperature (DragTempForcing.cpp)

    u*        = kappa * |u_h1| / ( ln(1.5 dz / z0) - psi_m )
    theta*    = theta1 * u*^2 / ( kappa * g * L )
    theta_s   = theta[k+1] - theta*/kappa * ( ln(1.5 dz / z0) - psi_h(1.5 dz / L) )
    theta_tgt = theta_s   + theta*/kappa * ( ln(0.5 dz / z0) - psi_h(0.5 dz / L) )
    bc_force  = - (theta_tgt - theta1) / (bc_forcing_time_factor * dt)
    src      -= bc_force * terrain_drag

- theta* follows from the definition L = u*^2 theta / (kappa g theta*), so the surface heat
  flux is set implicitly by the prescribed Obukhov length rather than by a flux input.
- Engineering choices that differ from the momentum side: u* is built from the drag cell
  itself rather than the neighbour, the cell above is always k+1, and the source is scaled by
  the signed marker. These are fine for terrain (marker is only ever +1 there) and only
  matter for overhangs / channel ceilings from ChannelBuilder.
- One item that looks like a slip rather than a choice: the heat-flux constants are read from
  `mo_gamma_m` / `mo_beta_m`, whereas `ABLWallFunction` reads `mo_gamma_h` / `mo_beta_h`.
  Harmless with defaults since all four default to the same values.

---

## 4. Code additions beyond the paper: stability

### 4.1 Atmospheric stability (Monin-Obukhov corrections)

Enabled by `ABL.wall_het_model = mol` with `ABL.monin_obukhov_length`, `ABL.kappa`,
`ABL.mo_beta_m`, `ABL.mo_gamma_m` (defaults 0.41, 16, 5). Functions in `MOData.cpp`:

    stable   (zeta > 0): psi_m = psi_h = -gamma * zeta
    unstable (zeta < 0): x = (1 - beta_m zeta)^(1/4)
                         psi_m = 2 ln((1+x)/2) + ln((1+x^2)/2) - 2 atan(x) + pi/2
                         psi_h = 2 ln((1 + sqrt(1 - beta_h zeta))/2)

These are evaluated once per level at zeta = 1.5 dz / L and 0.5 dz / L and enter the log
laws above. The paper's runs were neutral at the wall (no-slip, no stability functions);
the theory page's remark that "stability functions can be added in a straightforward
manner" is what these implement.

### 4.2 Numerical stability at start-up (the drag limiters)

The paper reports that increasing alpha_m beyond 1e3 caused numerical instabilities. The
code's defaults and limiters address the same problem for the explicit source term
du/dt = -c u with c = drag_coefficient/dz: forward Euler is monotone only for c dt <= 1.
At start-up the initial velocity is non-zero inside blanked cells, so a large c dt
overshoots and can blow up. Hence:

- `DragForcing.drag_coefficient` default 10 (not the paper's 1e3), and the docs recommend
  keeping the default "to avoid initial numerical instability".
- `DragForcing.use_original_drag_limiter` (default true) caps the velocity-normalised
  coefficient at `max_drag_coefficient` (default 1000) scaled by dz, and drops the 1/dz
  scaling for dz < 1 m (mirrors the paper's `max(1/Delta, 1)` for laboratory scales).
- `DragForcing.use_temporal_drag_limiter` (default false) caps the applied coefficient at
  1/dt so the relaxation can never exceed one time step (added in PR #1983 together with
  `max_drag_coefficient`, `is_laminar`, and ChannelBuilder options to initialise drag cells
  and skip velocity initialisation inside bodies).
- `DragForcing.bc_forcing_time_factor` (default 5) plays the same role for the wall-function
  relaxation: the log-law target is approached over 5 dt instead of 1 dt.
- `DragTempForcing` hard-codes its cap at 10/dz and has no temporal limiter.

### 4.3 Related terms with no counterpart in the paper

- Lateral sponge toward a 1-D RANS profile (`sponge_*`, quadratic ramp, `sponge_strength`).
  West/south distances are entered as negative numbers by convention.
- `terrain_damping`: sin^2 Rayleigh layers at the lateral boundaries (above
  `horizontal_abl_height`) and at the domain top, applied to w only, rate 1/`horizontal_tau`.

---

## 5. Convergence check: laminar immersed box (2026-09-05)

Case `test/test_files/terrain_box` (binary TerrainDrag + DragForcing, with
`use_temporal_drag_limiter = 1` so the explicit drag survives the 16 s cold-start step) versus
`test/test_files/immersed_terrain_box` (ImmersedTerrain + ImmersedDragForcing). Meshes
24x24x48, 48x48x96, 96x96x192 on the 1024 m cube, run to t = 90 s, inflow 1 m/s. Metric: mean
speed inside the box (true box extent x,y in [407, 593], z in [0, 200]).

| Region                     | Method  | dx = 42.7 | dx = 21.3 | dx = 10.7 | fitted order |
|----------------------------|---------|-----------|-----------|-----------|--------------|
| cells fully solid          | binary  | 1.20e-1   | 5.56e-2   | 2.05e-2   | 1.28 |
| cells fully solid          | partial | 1.18e-1   | 5.38e-2   | 2.59e-2   | 1.10 |
| interior window (430-570)  | binary  | 8.86e-2   | 3.25e-2   | 9.96e-3   | 1.58 |
| interior window (430-570)  | partial | 9.21e-2   | 3.16e-2   | 1.25e-2   | 1.44 |

Both methods are first order in the approach to zero inside the body; the partial fraction does
not improve it, as expected: in the laminar pathway the interior residual is pressure leakage through
the same relaxation term, and the fraction only alters the one-cell layer at the top of the box.
The gain from the partial fraction is in the surface layer (no staircase), which this metric does not
see. Without the temporal limiter the old explicit drag overshoots at the first step (C dt ~ 17) and
the CFL time step collapses to ~0.25 s; the exact-integration form in ImmersedDragForcing holds a
steady ~8 s step. Plot: `docs/notes/immersed_box_convergence_uniform.png`.

### 5.1 Why first order, and the implicit projection fix (2026-09-05)

Time-step and drag-coefficient sensitivity on the 48x48x96 mesh, t = 90 s, mean speed in fully
solid cells:

| Case                                   | steps | mean speed in solid |
|----------------------------------------|-------|---------------------|
| explicit drag, CFL dt (~9 s), Cd = 10  | 9     | 5.38e-2 |
| explicit drag, dt = 4.5, Cd = 10       | 20    | 3.20e-2 |
| explicit drag, dt = 2.25, Cd = 10      | 40    | 1.83e-2 |
| explicit drag, dt = 4.5, Cd = 100      | 20    | 3.15e-2 |
| implicit projection, CFL dt, Cd = 10   | 13    | 9.28e-3 |
| implicit projection, dt = 4.5          | 20    | 7.77e-3 |
| implicit projection, dt = 2.25         | 40    | 7.76e-3 |

With the explicit drag the residual halves when dt halves and ignores a tenfold larger Cd: the
drag zeroes the velocity, then the nodal projection re-injects dt grad(p)/rho inside the body,
because the projection treats the body as fluid. Under CFL dt ~ dx, hence first order in dx.

Fix (`ImmersedTerrain.implicit_projection = true`): the drag is applied implicitly with the
pressure. u^{n+1} = (u** - dt grad p / rho) / (1 + beta C dt), which makes the projection
coefficient sigma = dt / (rho (1 + beta C dt)); the MAC projection uses the same factor on faces.
Implemented in `incflo_apply_nodal_projection.cpp` and `icns_advection.cpp`, keyed on the
existence of the `terrain_drag_rate` field (= beta Cd / dz) that ImmersedTerrain declares.
ImmersedDragForcing then skips its explicit drag and keeps the wall model. Result: the interior
residual is independent of dt and six times smaller at the CFL step. The remaining residual is
grad(p)/(rho C) with C = Cd/dz, so it still scales with dz unless Cd is increased with
resolution (C ~ 1/dz^2 for second order), which the implicit form allows at no stability cost.

### 5.2 Convergence with the implicit projection (2026-09-05)

Same box case, meshes 24/48/96, t = 90 s. Mean speed in fully solid cells and in the interior
window; fitted order over the three levels.

| Series                                              | dx = 42.7 | dx = 21.3 | dx = 10.7 | order (solid) | order (interior) |
|-----------------------------------------------------|-----------|-----------|-----------|---------------|------------------|
| binary, explicit drag (with 1/dt limiter)           | 1.20e-1   | 5.56e-2   | 2.05e-2   | 1.28 | 1.58 |
| partial fraction, explicit drag                     | 1.18e-1   | 5.38e-2   | 2.59e-2   | 1.10 | 1.44 |
| partial fraction, implicit projection, Cd = 10      | 1.64e-2   | 9.28e-3   | 4.68e-3   | 0.91 | 1.02 |
| partial fraction, implicit projection, Cd = 10,20,40| 1.64e-2   | 4.12e-3   | 1.15e-3   | 1.92 | 2.03 |

Reading: the implicit projection removes the dt grad(p)/rho re-injection, cutting the residual by
5-10x at fixed Cd but leaving the spatial part grad(p)/(rho C) ~ dz. Doubling Cd with each
refinement (C ~ 1/dz^2) gives second order, which the implicit form allows at no stability cost.
Regression case: `test/test_files/immersed_terrain_box_implicit`.

### 5.3 Convergence with one AMR level on the surface (2026-09-05)

Same box, base meshes 24/48/96 with `amr.max_level = 1` and `FieldRefinement` on the terrain
surface band (`terrain_blank` grad_error 0.1 for the binary method, `terrain_mask` grad_error 0.5
for the new physics; the fine level covers 2.5-6.5 % of the domain). The x axis is the finest
spacing dx/2. t = 90 s, mean speed in fully solid cells and in the interior window.

| Series                                              | dx_f = 21.3 | dx_f = 10.7 | dx_f = 5.3 | order (solid) | order (interior) |
|-----------------------------------------------------|-------------|-------------|------------|---------------|------------------|
| binary, explicit drag (with 1/dt limiter)           | 5.55e-2     | 2.10e-2     | 8.30e-3    | 1.37 | 1.62 |
| partial fraction, explicit drag                     | 5.40e-2     | 2.66e-2     | 1.06e-2    | 1.18 | 1.51 |
| partial fraction, implicit projection, Cd = 10      | 9.16e-3     | 4.68e-3     | 2.57e-3    | 0.92 | 1.18 |
| partial fraction, implicit projection, Cd = 10,20,40| 9.16e-3     | 2.29e-3     | 6.82e-4    | 1.87 | 2.17 |

Consistency check: the AMR results at base 24 (finest 21.3 m) match the uniform 48 results of
Section 5.2 to within 2 % for every series (e.g. binary 5.55e-2 vs 5.56e-2, implicit 9.16e-3 vs
9.28e-3), so refining only the surface band reproduces the uniformly fine answer inside the body.
The orders are the same as on uniform meshes: first order for the explicit methods and for the
implicit projection at fixed Cd, second order for the implicit projection with Cd ~ 1/dz.
Regression cases: `test/test_files/immersed_terrain_box_amr` and `immersed_terrain_box_amr_implicit`
(tagging on `terrain_mask`; the `mask_terrain` derived sampling field is not used because it
requires the old `terrain_blank` int field). Plot: `docs/notes/immersed_box_convergence_amr.png`.

## 6. Diffusive flux at the terrain interface (2026-09-06)

The diffusion operator assembles face coefficients with `average_velocity_eta_to_faces` (six call
sites in `DiffusionOps.H` and `icns_diffusion.H`). At a fluid/solid face it therefore applies the
flux mu_eff (u_k - 0)/dx_f, since the interior is at rest. This is not a missing flux but a wall
stress with the wrong length scale: the no-slip flux for a wall at distance d1 is mu u_k/d1 (a
factor dx_f/d1, i.e. 2 for a wall at the face). With a turbulence model the wall model already
supplies u*^2, so the SGS flux across the interface counts the stress twice.

Option `ImmersedTerrain.interface_diffusion` (default `none`): ImmersedTerrain declares face fields
`terrain_diffusion_xf/yf/zf` and `diffusion::apply_immersed_interface` multiplies the face
coefficients by them at all six sites (so momentum, temperature, TKE and passive scalars are all
covered):

- `block`: factor min(1 - beta_L, 1 - beta_R); zero on faces touching the terrain. Turbulent
  pathway: the wall model alone carries stress and heat flux.
- `no_slip`: on fluid/solid faces (solid = fraction >= solid_threshold) factor dx_f/d1 with
  d1 = z_k - h on the bottom face (clamped to [0.1, 1] dz) and dx_f/2 elsewhere. Laminar pathway
  with finite viscosity: the discrete flux equals the no-slip flux at the true wall position.

Regression cases `immersed_terrain_box_viscous` and `immersed_terrain_box_viscous_amr` (nu = 1,
implicit projection, no_slip). Unit test `test_immersed_interface_diffusion.cpp` checks the factors
on the plateau: lateral wall 2, bottom face of the partial cell 32/12, blocked faces 0 and 0.875.

### 6.1 Viscous box, nu = 1 m^2/s (2026-09-06)

Laminar box with `transport.viscosity = 1.0` (cell Reynolds number U dx/nu = 5-43), implicit
projection, Cd = 10, t = 90 s, three interface treatments, uniform 24/48/96 and one AMR level.
Mean speed in fully solid cells (interior residual) and in the surface cells (first fluid layer,
`terrain_mask` = 2). No exact reference exists for this flow; the table shows trends only.

| Grid | Metric        | Mode    | coarse   | medium   | fine     | fitted slope |
|------|---------------|---------|----------|----------|----------|--------------|
| uni  | solid cells   | none    | 1.64e-2  | 7.89e-3  | 4.11e-3  | 1.00 |
| uni  | solid cells   | no_slip | 1.66e-2  | 7.99e-3  | 4.31e-3  | 0.97 |
| uni  | solid cells   | block   | 1.66e-2  | 8.04e-3  | 4.28e-3  | 0.98 |
| uni  | surface cells | none    | 0.728    | 0.631    | 0.443    | 0.36 |
| uni  | surface cells | no_slip | 0.720    | 0.587    | 0.361    | 0.50 |
| uni  | surface cells | block   | 0.737    | 0.693    | 0.589    | 0.16 |
| amr  | solid cells   | none    | 7.90e-3  | 4.11e-3  | 2.19e-3  | 0.93 |
| amr  | solid cells   | no_slip | 7.99e-3  | 4.31e-3  | 2.24e-3  | 0.92 |
| amr  | solid cells   | block   | 8.04e-3  | 4.28e-3  | 2.25e-3  | 0.92 |
| amr  | surface cells | none    | 0.631    | 0.443    | 0.292    | 0.56 |
| amr  | surface cells | no_slip | 0.587    | 0.361    | 0.205    | 0.76 |
| amr  | surface cells | block   | 0.693    | 0.589    | 0.607    | 0.10 |

(uniform dx = 42.7 / 21.3 / 10.7 m; AMR finest dx = 21.3 / 10.7 / 5.3 m)

Reading: the interior residual is unaffected by the interface treatment, as expected. The near-wall
speed converges toward the resolved no-slip layer fastest with `no_slip` (slope 0.76 on the AMR
series, heading to the first-order behavior u(d1) ~ tau_w d1/mu of a resolved layer), slower with
`none` (wall stress too weak by dx/d1), and not at all with `block` (no viscous wall stress, a slip
wall in the laminar pathway; intended for the turbulent pathway where the wall model supplies the
stress). AMR at base 24 again matches uniform 48 for every entry. Plot:
`docs/notes/immersed_box_viscous_interface.png`.

## 7. Wall-model consistency: manufactured log-law test (2026-09-08)

`unit_tests/wind_energy/test_immersed_wall_mms.cpp`. An exact neutral log-law velocity,
tangential to the terrain and a function of the true closest-point distance, is imposed above
(a) a plane slope of 20 degrees and (b) a Gaussian ridge (A = 200 m, sigma = 200 m, max slope 31
degrees). The wall-model geometry of `ImmersedWallModel.H` recovers u* and the log-law target
in every surface cell; the table gives the mean relative error over the surface cells.

| Terrain | Method                            | n = 32   | n = 64   | n = 128  | behaviour |
|---------|-----------------------------------|----------|----------|----------|-----------|
| plane   | cell_offset                       | 0.340    | 0.365    | 0.364    | O(1), staircase + slope |
| plane   | terrain_height                    | 0.331    | 0.356    | 0.355    | O(1), vertical vs normal distance |
| plane   | surface_normal, nominal d2        | 2.4e-2   | 3.1e-2   | 3.5e-2   | O(1), reference distance inconsistent |
| plane   | surface_normal, actual d2         | 1.5e-16  | 2.5e-16  | 5.3e-16  | exact |
| ridge   | cell_offset                       | 0.366    | 0.368    | 0.354    | O(1) |
| ridge   | terrain_height                    | 0.356    | 0.354    | 0.341    | O(1) |
| ridge   | surface_normal, nominal d2        | 2.5e-2   | 2.9e-2   | 3.2e-2   | O(1) |
| ridge   | surface_normal, actual d2         | 6.2e-4   | 3.4e-4   | 2.1e-4   | converges, order 0.77 (32->128) |
| ridge   | surface_normal, actual d2, center | 7.6e-4   | 4.0e-4   | 2.4e-4   | converges, order 0.83 |

(relative error in u*; target errors follow the same pattern, e.g. plane actual d2 1e-14.)

Reading: the two new options are `ImmersedDragForcing.reference_distance = actual` (d2 from the
sampled reference cell's own normal distance) and `ImmersedTerrain.drag_weight = center`. With
`actual` the surface_normal wall model is exact on a plane and converges on a curved surface; the
remaining per-cell error is the tangent-plane approximation of the distance on a curved surface,
~ dx / ln(dx/z0), hence an order a little below 1 per cell (second order in the volume average).
The face-based methods and the nominal normal method carry a slope error that does not decrease
with the mesh (33-36 % and 2-3 % of u* respectively at 20-31 degrees). The debug dump showed the
ridge error spread smoothly over all cells (none above 1e-3), i.e. no outliers from the
nearest-cell reference rounding.
