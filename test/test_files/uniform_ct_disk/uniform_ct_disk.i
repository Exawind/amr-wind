time.stop_time               =   -100.0     # Max (simulated) time to evolve
time.max_step                =   10         # Max number of time steps

time.fixed_dt         =   0.1      # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor

io.KE_int = 0
io.outputs = actuator_src_term
time.plot_interval            =  10       # Steps between plot files
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
Actuator.type = UniformCtDisk

Actuator.UniformCtDisk.rotor_diameter = 126.0
Actuator.UniformCtDisk.disk_normal = 10 0.0 0.0
Actuator.UniformCtDisk.disk_center = 0.0 0.0 0.0
Actuator.UniformCtDisk.num_force_points = 5
Actuator.UniformCtDisk.thrust_coeff = 0.7
Actuator.UniformCtDisk.epsilon = 10.0
Actuator.UniformCtDisk.density = 1.225
Actuator.UniformCtDisk.diameters_to_sample = 1.0
Actuator.UniformCtDisk.num_vel_points_r = 3
Actuator.UniformCtDisk.num_vel_points_t = 3

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
nodal_proj.verbose = 0

nodal_proj.mg_rtol = 1.0e-6
nodal_proj.mg_atol = 1.0e-12
mac_proj.mg_rtol = 1.0e-6
mac_proj.mg_atol = 1.0e-12

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              Sampling                 #
#.......................................#
incflo.post_processing = sampling line_sampling averaging
# Frequency of output for the data
sampling.output_frequency = 1000
sampling.labels = p_hub
# Fields to output
sampling.fields = velocity velocity_mean velocity_reynolds_stress
# Definitions for each probe
sampling.p_hub.type = PlaneSampler
sampling.p_hub.axis2 = 0.0 0.0 630.0
sampling.p_hub.axis1 = 630 0.0 0.0
sampling.p_hub.origin = -315.0 0.0 -315.0
sampling.p_hub.num_points = 128 128
sampling.p_hub.normal = 0.0 1.0 0.0

line_sampling.output_frequency = 10
line_sampling.labels = line0 line1
line_sampling.fields = velocity velocity_mean velocity_reynolds_stress

line_sampling.line0.type       = LineSampler
line_sampling.line0.num_points = 160
line_sampling.line0.start      = 0.0 0.0 -126.0
line_sampling.line0.end        = 0.0 0.0  126.0

line_sampling.line1.type       = LineSampler
line_sampling.line1.num_points = 160
line_sampling.line1.start      = 126.0 0.0 -126.0
line_sampling.line1.end        = 126.0 0.0  126.0

# The time averaging
averaging.type = TimeAveraging
averaging.labels = means  stress

averaging.averaging_window = 10.0
averaging.averaging_start_time = 0.0

averaging.means.fields = velocity
averaging.means.averaging_type = ReAveraging

averaging.stress.fields = velocity
averaging.stress.averaging_type = ReynoldsStress
