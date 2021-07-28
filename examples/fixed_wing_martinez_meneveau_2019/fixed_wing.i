#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   5.0     # Max (simulated) time to evolve
time.max_step                =  -500          # Max number of time steps

tagging.labels = static
tagging.static.type = CartBoxRefinement
tagging.static.static_refinement_def = static_box.txt

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.005        # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
io.KE_int = 0
io.outputs = actuator_src_term
time.plot_interval            =  100       # Steps between plot files
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

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          FIXED WING INPUTS            #
#.......................................#
incflo.physics = FreeStream Actuator
Actuator.labels = F1
Actuator.type = FixedWingLine
Actuator.FixedWingLine.num_points = 600
Actuator.FixedWingLine.epsilon_chord = 1.0 1.0 1.0
Actuator.FixedWingLine.epsilon = .001 .001 .001
Actuator.FixedWingLine.pitch = 6.0
Actuator.FixedWingLine.span_locs = 0.0 1.0
Actuator.FixedWingLine.chord = .08 .08
Actuator.FixedWingLine.airfoil_table = NACA64_A17.dat
Actuator.FixedWingLine.airfoil_type = openfast
Actuator.F1.start = 0.0 -0.5 0.0
Actuator.F1.end =   0.0  0.5 0.0
Actuator.F1.output_frequency = 10

ICNS.source_terms = ActuatorForcing

amr.n_cell              = 96 64 64    # Grid cells at coarsest AMRlevel
amr.max_level           = 1           # Max AMR level in hierarchy 
geometry.prob_lo        = -2. -2. -2.
geometry.prob_hi        =  4.  2.  2.

geometry.is_periodic    =   0   0   0   # Periodicity x y z (0/1)

# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.0
xlo.velocity = 1.0 0.0 0.0
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
