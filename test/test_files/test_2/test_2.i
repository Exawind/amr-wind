#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
stop_time               =   -1          # Max (simulated) time to evolve
max_step                =   20        # Max number of time steps
steady_state            =   0           # Steady-state solver? 
incflo.steady_state_tol =   1.e-5       # Tolerance for steady-state

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
#incflo.fixed_dt         =   0.01        # Use this constant dt if > 0
incflo.cfl              =   0.45         # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
#amr.KE_int = 1
amr.plot_int            =   20       # Steps between plot files
amr.plot_per            =   -1          # Steps between plot files
amr.check_int           =  -100       # Steps between checkpoint files
amr.restart             =   ""          # Checkpoint to restart from 
amr.plt_tracer = 1
amr.plt_eta = 1
amr.plt_p = 1
amr.plt_divu = 1
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity          =   0.  0. -9.80665  # Gravitational force (3D)
incflo.ro_0             = 1.225          # Reference density 

incflo.fluid_model      =   "newtonian" # Fluid model (rheology)
incflo.mu               =   1.849e-5      # Dynamic viscosity coefficient

abl.use_boussinesq = 1 
abl.coriolis_effect = 1
abl.abl_forcing = 1
abl.sgs_model = 1

incflo.advect_tracer = 1
incflo.ntrac = 1 

incflo.probtype = 35
incflo.ic_u = 6.76
incflo.ic_v = -5.94

abl.ntemperature = 3
abl.temperature_heights = 650.0 750.0 1000.0
abl.temperature_values = 300.0 308.0 308.75

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              =  32 32 32   # Grid cells at coarsest AMRlevel
amr.max_level           =   0           # Max AMR level in hierarchy 

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.       0.     0.  # Lo corner coordinates
geometry.prob_hi        =   1000.  1000.  1000.  # Hi corner coordinates
geometry.is_periodic    =   1   1   0   # Periodicity x y z (0/1)

incflo.delp             =   0.  0.  0.  # Prescribed (cyclic) pressure gradient

# Boundary conditions
zlo.type =   "wall_model"
zlo.tracer = 0.0

zhi.type =   "slip"
zhi.tracer = 0.003 # tracer is used to specify potential temperature gradient

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              VERBOSITY                #
#.......................................#
incflo.verbose          =   2           # incflo_level
diffusion.verbose       =   0           # DiffusionEquation
mac.verbose             =   0           # MacProjector

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              DEBUGGING                #
#.......................................#
amrex.fpe_trap_invalid  =   0           # Trap NaNs


diffusion.verbose = 1
diffusion.mg_verbose = 1
diffusion.mg_cg_verbose = 1
diffusion.mg_rtol = 1.0e-8
diffusion.mg_atol = 1.0e-8

