#include <type_traits>
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   300.00     # Max (simulated) time to evolve
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
BodyForce.magnitude = 0.0 6e-2 0.0
incflo.physics = ChannelFlow
ChannelFlow.density = 1.0
ChannelFlow.Laminar = true
ChannelFlow.flow_direction = 1
ChannelFlow.normal_direction = 2
ChannelFlow.Turbulent_DNS = false
ChannelFlow.Mean_Velocity = 1.0

io.output_default_variables = 1

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 16 32 32    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.0  0.0  0.0  # Lo corner coordinates
geometry.prob_hi        =   1.0  6.0  1.0  # Hi corner coordinates
geometry.is_periodic    =   1    1    0   # Periodicity x y z (0/1)

geometry.mesh_mapping = ChannelFlowMap
ChannelFlowMap.beta = 0 0 3

velocity_diffusion.use_tensor_operator = false
incflo.diffusion_type   = 2
    
# Boundary conditions
zlo.type =   "no_slip_wall"
zhi.type =   "no_slip_wall"

incflo.verbose  = 0

diffusion.max_coarsening_level = 0

nodal_proj.mg_atol = 1.0e-09
nodal_proj.verbose = 0
nodal_proj.bottom_solver = hypre
nodal_proj.max_coarsening_level = 0

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
hypre.recompute_preconditioner = 0
hypre.verbose = 0
hypre.bamg_verbose = 0
hypre.num_krylov = 20
hypre.bamg_max_levels = 4
hypre.bamg_num_sweeps = 1
#
