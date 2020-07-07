time.max_step                = 1000000
time.stop_time               = 5.0           # Max (simulated) time to evolve

incflo.use_godunov                     = true

time.fixed_dt = 0.05
time.cfl              = 0.49           # CFL factor
time.init_shrink      = 1.0

time.plot_interval            =   140         # Steps between plot files
time.checkpoint_interval           =  -100         # Steps between checkpoint files

transport.viscosity = 1.0
transport.laminar_prandtl = 0.33333333333333337
turbulence.model = Laminar

amr.max_level           =   0
amr.n_cell              =   16 16 16     # Grid cells at coarsest AMRlevel
amr.max_grid_size       =   16 16 16  # Max grid size at AMR levels
amr.blocking_factor     =   8            # Blocking factor for grids

geometry.prob_lo        = -3.14159265358979323 -3.14159265358979323 -3.14159265358979323 # Lo corner coordinates
geometry.prob_hi        =  3.14159265358979323  3.14159265358979323  3.14159265358979323 # Hi corner coordinates
  
geometry.is_periodic    =   1   1   1   # Periodicity x y z (0/1)

incflo.probtype         =  0
incflo.gravity          = 0. 0. 0.

incflo.physics = MMS
MMS.masa_name = "navierstokes_3d_incompressible_homogeneous"
ICNS.source_terms = MMSForcing

incflo.diffusion_type   = 2             # 0 = Explicit, 1 = Crank-Nicolson, 2 = Implicit

incflo.verbose          =   1           # incflo_level
mac_proj.verbose        =   0           # MAC Projector
nodal_proj.verbose      =   0           # Nodal Projector
diffusion.mg_verbose    =   0           # Diffusion

mac_proj.mg_rtol        = 1.e-12
