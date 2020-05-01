time.max_step                = 20

time.fixed_dt         = 0.1
time.cfl              = 0.9           # CFL factor
time.init_shrink      = 1.0

time.plot_interval            =   10          # Steps between plot files
time.checkpoint_interval           =  -100         # Steps between checkpoint files

transport.viscosity = 0.001
transport.laminar_prandtl = 1.0
transport.turbulent_prandtl = 1.0
turbulence.model = Laminar

amr.max_level           =   2
time.regrid_interval          =   2

amr.n_cell              =   16 16 32    # Grid cells at coarsest AMRlevel
amr.max_grid_size       =   8 8 8       # Max grid size at AMR levels
amr.blocking_factor     =   8           # Blocking factor for grids

geometry.prob_lo        =  0.  0.  0.   # Lo corner coordinates
geometry.prob_hi        =  0.5 0.5 1.0  # Hi corner coordinates

geometry.is_periodic    =   1   1   0   # Periodicity x y z (0/1)

# Boundary
zlo.type                = "sw"
zhi.type                = "sw"

incflo.probtype         =  0
incflo.physics = RayleighTaylor
ICNS.source_terms = DensityBuoyancy
RayleighTaylor.rho_lo = 0.5
RayleighTaylor.rho_hi = 2.0

incflo.gravity          = 0. 0. -0.3

incflo.gradrhoerr       = 0.1

amr.plt_ccse_regtest    =  1
amr.plt_vort            =  1

incflo.use_godunov      = true
incflo.constant_density = false
incflo.advect_tracer    = true

incflo.diffusion_type   = 1             # 0 = Explicit, 1 = Crank-Nicolson, 2 = Implicit

incflo.verbose          =   0           # incflo_level

mac_proj.mg_rtol        = 1.e-12
