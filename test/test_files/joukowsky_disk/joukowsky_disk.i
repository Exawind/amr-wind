time.stop_time               =   -100.0     # Max (simulated) time to evolve
time.max_step                =   10         # Max number of time steps

time.fixed_dt         =   0.1      # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor

io.KE_int = 0
io.outputs = actuator_src_term
time.plot_interval            =  1       # Steps between plot files
time.checkpoint_interval      =  -1       # Steps between checkpoint files

ConstValue.density.value = 1.0
ConstValue.velocity.value = 8.0 0.0 0.0

incflo.use_godunov = 1
incflo.godunov_type = "ppm_nolim"
incflo.do_initial_proj = 1
incflo.initial_iterations = 3
transport.viscosity = 1.0e-5
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model               = Smagorinsky
Smagorinsky_coeffs.Cs          = 0.16

incflo.physics = FreeStream Actuator
Actuator.labels = WTG01
Actuator.type = JoukowskyDisk

Actuator.JoukowskyDisk.rotor_diameter = 126.0
Actuator.JoukowskyDisk.disk_center = 0.0 0.0 0.0
Actuator.JoukowskyDisk.num_blades = 3
Actuator.JoukowskyDisk.yaw = 270.0 # degrees (yaw is relative to north which defaults to {0,1,0})
Actuator.JoukowskyDisk.sample_yaw = 270.0 # set velocity sampling to be in the normal flow direction
Actuator.JoukowskyDisk.thrust_coeff = 0.0 0.7 1.2
Actuator.JoukowskyDisk.wind_speed = 0.0 10.0 12.0
Actuator.JoukowskyDisk.angular_velocity = 0.0 1.0 1.5
Actuator.JoukowskyDisk.epsilon = 20.0
Actuator.JoukowskyDisk.density = 1.225
Actuator.JoukowskyDisk.diameters_to_sample = 1.0
Actuator.JoukowskyDisk.num_points_r = 5
Actuator.JoukowskyDisk.num_points_t = 3

ICNS.source_terms = ActuatorForcing

amr.n_cell              = 64 64 64   # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy
geometry.prob_lo        =   -315.0 -315.0 -315.0
geometry.prob_hi        =   315.0  315.0  315.0

geometry.is_periodic    =   0   0   0   # Periodicity x y z (0/1)

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        Mesh refinement                #
#.......................................#
tagging.labels = refine0
tagging.refine0.type = GeometryRefinement
tagging.refine0.shapes = c0 c1

tagging.refine0.c0.type = cylinder
tagging.refine0.c0.start = -256.0 0 0
tagging.refine0.c0.end = 2600.0 0 0
tagging.refine0.c0.outer_radius = 140.0 

tagging.refine0.c1.type = cylinder
tagging.refine0.c1.start = -189.0 0 0
tagging.refine0.c1.end = 2646.0 0 0
tagging.refine0.c1.outer_radius = 70.0

# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.0
xlo.velocity = 8.0 0.0 0.0
xhi.type = "pressure_outflow"
ylo.type =   "slip_wall"
yhi.type =   "slip_wall"
zlo.type =   "slip_wall"
zhi.type =   "slip_wall"

incflo.verbose          =   0          # incflo_level
