#include <type_traits>
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   100.00     # Max (simulated) time to evolve
time.max_step                =   -1000       # Max number of time steps
incflo.initial_iterations    = 0
incflo.do_initial_proj       = 0

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.03       # Use this constant dt if > 0
time.cfl              =   0.5       # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval            =  1000       # Steps between plot files
time.checkpoint_interval      =  -100       # Steps between checkpoint files

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.density        =  1.0             # Reference density
incflo.use_godunov = 0
transport.viscosity = 0.005
turbulence.model = Laminar

ICNS.source_terms = BodyForce
BodyForce.magnitude = 6e-2 0 0
incflo.physics = ChannelFlow
ChannelFlow.density = 1.0
ChannelFlow.Laminar = true
ChannelFlow.Turbulent_DNS = false
ChannelFlow.Mean_Velocity = 1.0

io.output_default_variables = 1

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 32 32 16    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.0  0.0  0.0  # Lo corner coordinates
geometry.prob_hi        =   6.0  1.0  1.0  # Hi corner coordinates
geometry.is_periodic    =   1   0   1   # Periodicity x y z (0/1)

geometry.mesh_mapping   = ChannelFlowScaling
ChannelFlowScaling.beta = 0 3 0
ChannelFlowScaling.do_map = 0 1 0

velocity_diffusion.use_tensor_operator = false
incflo.diffusion_type   = 2
    
# Boundary conditions
ylo.type =   "no_slip_wall"
yhi.type =   "no_slip_wall"

diffusion.bottom_solver = hypre
diffusion.max_coarsening_level = 0

incflo.verbose  = 0
nodal_proj.mg_atol = 1.0e-09
nodal_proj.verbose = 0
nodal_proj.bottom_solver = hypre
nodal_proj.max_coarsening_level = 0
#
mac_proj.mg_rtol = 1.0e-11
mac_proj.mg_atol = 1.0e-09
mac_proj.do_semicoarsening = true
mac_proj.bottom_solver = hypre
mac_proj.bottom_verbose = 0
mac_proj.max_coarsening_level = 0
mac_proj.bottom_rtol = 1.0e-12
mac_proj.bottom_atol = 1.0e-16
#
hypre.hypre_solver = GMRES
hypre.hypre_preconditioner = BoomerAMG
hypre.verbose = 0
hypre.bamg_verbose = 0
hypre.num_krylov = 20
hypre.bamg_max_levels = 4
hypre.bamg_num_sweeps = 1
#
