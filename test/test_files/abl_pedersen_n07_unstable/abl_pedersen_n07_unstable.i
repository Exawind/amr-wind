#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =  82800.0     # Max (simulated) time to evolve
time.max_step                =   -1          # Max number of time steps
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.5        # Use this constant dt if > 0
time.cfl              =   0.95       # CFL factor
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval            =  1000       # Steps between plot files
time.checkpoint_interval      =  10000       # Steps between checkpoint files
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity        =  0.0  0.0 -9.81  # Gravitational force (3D)
incflo.density        =  1.0             # Reference density
incflo.use_godunov = 1
incflo.use_limiter = 0
transport.viscosity = 0.0
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = OneEqKsgsM84
incflo.physics = ABL
ICNS.source_terms = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing
TKE.source_terms = KsgsM84Src
BoussinesqBuoyancy.reference_temperature = 290.0
CoriolisForcing.east_vector = 1.0 0.0 0.0
CoriolisForcing.north_vector = 0.0 1.0 0.0
CoriolisForcing.latitude = 90.0
CoriolisForcing.rotational_time_period = 125663.706143592
GeostrophicForcing.geostrophic_wind = 10.0 0.0 0.0
incflo.velocity = 10.0 0.0 0.0
ABL.reference_temperature = 290.0
ABL.temperature_heights = 0.0 2000.0
ABL.temperature_values = 290.0 296.0
ABL.perturb_temperature = true
ABL.cutoff_height = 50.0
ABL.perturb_velocity = true
ABL.perturb_ref_height = 50.0
ABL.Uperiods = 4.0
ABL.Vperiods = 4.0
ABL.deltaU = 1.0
ABL.deltaV = 1.0
ABL.kappa = .41
ABL.surface_roughness_z0 = 0.01
ABL.surface_temp_flux = 0.005
ABL.normal_direction = 2
ABL.stats_output_frequency = 1
#ABL.stats_output_format = netcdf


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 36 36 24    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy
amr.blocking_factor     = 4
tagging.static_refinement = false
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.       0.     0.  # Lo corner coordinates
geometry.prob_hi        =   3000.  3000.  2000.  # Hi corner coordinates
geometry.is_periodic    =   1   1   0   # Periodicity x y z (0/1)
# Boundary conditions
zlo.type =   "wall_model"
zlo.temperature_type = "wall_model"
zhi.type =   "slip_wall"
zhi.temperature_type = "fixed_gradient"
zhi.temperature = 0.003
zlo.tke_type = "fixed_gradient"
incflo.verbose          =   0
# MLMG settings
nodal_proj.mg_rtol = 1.0e-6
nodal_proj.mg_atol = 1.0e-12
mac_proj.mg_rtol = 1.0e-6
mac_proj.mg_atol = 1.0e-12
diffusion.mg_rtol = 1.0e-6
diffusion.mg_atol = 1.0e-12
temperature_diffusion.mg_rtol = 1.0e-10
temperature_diffusion.mg_atol = 1.0e-13

# Activate data probe sampling
incflo.post_processing = sampling

# Frequency of output for the data
sampling.output_frequency = 5
# Type of probes to output
sampling.labels = l_v11 l_v12 l_v13 l_v21 l_v22 l_v23 l_v31 l_v32 l_v33 p_h p_v1 p_v2

# Fields to output
sampling.fields = velocity temperature

# Definitions for each probe
sampling/l_v11.type = LineSampler
sampling/l_v11.num_points = 192
sampling/l_v11.start = 500.0 500.0 5.208333333333333
sampling/l_v11.end = 500.0 500.0 1994.7916666666667

sampling/l_v12.type = LineSampler
sampling/l_v12.num_points = 192
sampling/l_v12.start = 500.0 1500.0 5.208333333333333
sampling/l_v12.end = 500.0 1500.0 1994.7916666666667

sampling/l_v13.type = LineSampler
sampling/l_v13.num_points = 192
sampling/l_v13.start = 500.0 2500.0 5.208333333333333
sampling/l_v13.end = 500.0 2500.0 1994.7916666666667

sampling/l_v21.type = LineSampler
sampling/l_v21.num_points = 192
sampling/l_v21.start = 1500.0 500.0 5.208333333333333
sampling/l_v21.end = 1500.0 500.0 1994.7916666666667

sampling/l_v22.type = LineSampler
sampling/l_v22.num_points = 192
sampling/l_v22.start = 1500.0 1500.0 5.208333333333333
sampling/l_v22.end = 1500.0 1500.0 1994.7916666666667

sampling/l_v23.type = LineSampler
sampling/l_v23.num_points = 192
sampling/l_v23.start = 1500.0 2500.0 5.208333333333333
sampling/l_v23.end = 1500.0 2500.0 1994.7916666666667

sampling/l_v31.type = LineSampler
sampling/l_v31.num_points = 192
sampling/l_v31.start = 2500.0 500.0 5.208333333333333
sampling/l_v31.end = 2500.0 500.0 1994.7916666666667

sampling/l_v32.type = LineSampler
sampling/l_v32.num_points = 192
sampling/l_v32.start = 2500.0 1500.0 5.208333333333333
sampling/l_v32.end = 2500.0 1500.0 1994.7916666666667

sampling/l_v33.type = LineSampler
sampling/l_v33.num_points = 192
sampling/l_v33.start = 2500.0 2500.0 5.208333333333333
sampling/l_v33.end = 2500.0 2500.0 1994.7916666666667

sampling/p_h.type = PlaneSampler
sampling/p_h.axis1 = 1.0 0.0 0.0
sampling/p_h.axis2 = 0.0 1.0 0.0
sampling/p_h.origin = 0.0 0.0 50.0
sampling/p_h.num_points = 288 288
sampling/p_h.normal = 0.0 0.0 1.0
sampling/p_h.offsets = 0.0 50.0 150.0

sampling/p_v1.type = PlaneSampler
sampling/p_v1.axis1 = 0.0 0.0 1.0
sampling/p_v1.axis2 = 1.0 0.0 0.0
sampling/p_v1.origin = 0.0 1500.0 0.0
sampling/p_v1.num_points = 192 288

sampling/p_v2.type = PlaneSampler
sampling/p_v2.axis1 = 0.0 0.0 1.0
sampling/p_v2.axis2 = 0.0 -1.0 0.0
sampling/p_v2.origin = 1500.0 3000.0 0.0
sampling/p_v2.num_points = 192 288
