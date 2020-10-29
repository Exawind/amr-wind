#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   6     # Max (simulated) time to evolve
time.max_step                =   1000          # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.001        # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval            =  10       # Steps between plot files
time.checkpoint_interval      =  -100       # Steps between checkpoint files

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#

incflo.use_godunov = 1
incflo.godunov_type = "weno"
transport.viscosity = 1.
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Laminar 

incflo.physics = MultiPhase DamBreak 
MultiPhase.density_fluid1=100
MultiPhase.density_fluid2=1
ICNS.source_terms = DensityBuoyancy

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 64 16 64    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy 
amr.blocking_factor     = 8

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0   0.   0.  # Lo corner coordinates
geometry.prob_hi        =   0.5   0.125  0.5  # Hi corner coordinates
geometry.is_periodic    =   0   1   0   # Periodicity x y z (0/1)

xlo.type =   "slip_wall"
xhi.type =   "slip_wall"
zlo.type =   "slip_wall"
zhi.type =   "slip_wall"

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              VERBOSITY                #
#.......................................#
incflo.verbose          =   0          # incflo_level

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              DEBUGGING                #
#.......................................#
incflo.verbose=1
nodal_proj.mg_rtol= 1.0e-6
nodal_proj.mg_atol= 1.0e-10
diffusion.mg_rtol = 1.0e-6
diffusion.mg_atol = 1.0e-10
mac_proj.mg_rtol = 1.0e-6
mac_proj.mg_atol = 1.0e-10
