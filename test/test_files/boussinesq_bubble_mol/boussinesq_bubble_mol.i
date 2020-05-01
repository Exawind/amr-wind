time.max_step                = 140 
time.stop_time               = 100.0           # Max (simulated) time to evolve

incflo.use_godunov      = false

time.fixed_dt = 0.2
time.cfl              = 0.49           # CFL factor
time.init_shrink      = 1.0

time.plot_interval            =   140         # Steps between plot files
time.checkpoint_interval           =  -100         # Steps between checkpoint files

transport.viscosity = 0.00001
transport.laminar_prandtl = 0.33333333333333337
turbulence.model = Laminar

amr.max_level           =   1
amr.n_cell              =   32 32 64     # Grid cells at coarsest AMRlevel
amr.max_grid_size       =   16 16 16  # Max grid size at AMR levels
amr.blocking_factor     =   8            # Blocking factor for grids

geometry.prob_lo        =  0.  0.  0.   # Lo corner coordinates
geometry.prob_hi        =  0.5 0.5 1.0  # Hi corner coordinates

geometry.is_periodic    =   1   1   0   # Periodicity x y z (0/1)

incflo.gravity          = 0. 0. 0.5

incflo.physics = BoussinesqBubble
ICNS.source_terms = BoussinesqBuoyancy
BoussinesqBuoyancy.thermal_expansion_coeff = 1.0
BoussinesqBuoyancy.reference_temperature = 0.0

zlo.type                = "sw"
zhi.type                = "sw"

incflo.gradrhoerr       = 0.1

amr.plt_ccse_regtest    =  1
amr.plt_vort            =  1

incflo.diffusion_type   = 2             # 0 = Explicit, 1 = Crank-Nicolson, 2 = Implicit

incflo.verbose          =   1           # incflo_level
mac_proj.verbose        =   0           # MAC Projector
nodal_proj.verbose      =   0           # Nodal Projector
diffusion.mg_verbose    =   0           # Diffusion

mac_proj.mg_rtol        = 1.e-12
