#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   -100.0     # Max (simulated) time to evolve
time.max_step                =   10          # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.05        # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
io.KE_int = 0
io.outputs = actuator_src_term
time.plot_interval            =  50       # Steps between plot files
time.checkpoint_interval      =  -1       # Steps between checkpoint files

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
ConstValue.density.value = 1.0
ConstValue.velocity.value = 10.0 0.0 0.0

incflo.use_godunov = 1
incflo.do_initial_proj = 1
incflo.initial_iterations = 3
transport.viscosity = 1.0e-5
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Laminar

incflo.physics = FreeStream Actuator
Actuator.labels = F1
Actuator.type = FlatPlateLine
Actuator.FlatPlateLine.num_points = 21
Actuator.FlatPlateLine.epsilon = 3.0 3.0 3.0
Actuator.FlatPlateLine.pitch = 4.0
Actuator.F1.start = 0.0 -4.0 0.0
Actuator.F1.end = 0.0 4.0 0.0
Actuator.F1.output_frequency = 10

ICNS.source_terms = ActuatorForcing


amr.n_cell              = 64 32 32    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy 
geometry.prob_lo        =   -16.0 -16.0 -16.0
geometry.prob_hi        =   48.0 16.0 16.0
geometry.is_periodic    =   0   0   0   # Periodicity x y z (0/1)

# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.0
xlo.velocity = 10.0 0.0 0.0
xhi.type = "pressure_outflow"
ylo.type =   "slip_wall"
yhi.type =   "slip_wall"
zlo.type =   "slip_wall"
zhi.type =   "slip_wall"

incflo.verbose          =   0          # incflo_level
nodal_proj.verbose = 0

nodal_proj.mg_rtol = 1.0e-6
nodal_proj.mg_atol = 1.0e-12
mac_proj.mg_rtol = 1.0e-6
mac_proj.mg_atol = 1.0e-12
