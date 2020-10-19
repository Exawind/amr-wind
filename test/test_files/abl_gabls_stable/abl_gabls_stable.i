#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               = 28800.0     # Max (simulated) time to evolve
time.max_step                =   -1          # Max number of time steps
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.25        # Use this constant dt if > 0
time.cfl              =   0.95       # CFL factor
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
io.KE_int = -1
time.plot_interval            =  1000       # Steps between plot files
time.checkpoint_interval      =  10000       # Steps between checkpoint files
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity        =  0.0  0.0 -9.81  # Gravitational force (3D)
incflo.density        =  1.3223            # Reference density
incflo.use_godunov = 1
incflo.use_limiter = 0
transport.viscosity = 0.0
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = OneEqKsgsM84
incflo.physics = ABL
ICNS.source_terms = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing
TKE.source_terms = KsgsM84Src
BoussinesqBuoyancy.reference_temperature = 263.5
CoriolisForcing.east_vector = 1.0 0.0 0.0
CoriolisForcing.north_vector = 0.0 1.0 0.0
CoriolisForcing.latitude = 90.0
CoriolisForcing.rotational_time_period = 90405.5439881955
GeostrophicForcing.geostrophic_wind = 8.0 0.0 0.0
incflo.velocity = 8.0 0.0 0.0
ABL.reference_temperature = 263.5
ABL.temperature_heights = 0.0 100 400.0
ABL.temperature_values = 265.0 265.0 268.0
ABL.perturb_temperature = true
ABL.cutoff_height = 50.0
ABL.perturb_velocity = true
ABL.perturb_ref_height = 50.0
ABL.Uperiods = 4.0
ABL.Vperiods = 4.0
ABL.deltaU = 1e-5
ABL.deltaV = 1e-5
ABL.kappa = .40
ABL.mo_gamma_m = 4.8
ABL.mo_gamma_h = 7.8
ABL.surface_roughness_z0 = 0.1
ABL.surface_temp_rate = -0.25
ABL.surface_temp_init = 265.0
ABL.normal_direction = 2
ABL.stats_output_frequency = 2
#ABL.stats_output_format = netcdf
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 48 48 48 # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy
tagging.static_refinement = false
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.       0.     0.  # Lo corner coordinates
geometry.prob_hi        =   400.   400.   400.  # Hi corner coordinates
geometry.is_periodic    =   1   1   0   # Periodicity x y z (0/1)
# Boundary conditions
zlo.type =   "wall_model"
zlo.temperature_type = "wall_model"
zhi.type =   "slip_wall"
zhi.temperature_type = "fixed_gradient"
zhi.temperature = 0.01
zlo.tke_type = "zero_gradient"
#zhi.tke_type = "zero_gradient"
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
sampling.output_frequency = 2
# Type of probes to output
sampling.labels = l_v11 l_v12 l_v13 l_v21 l_v22 l_v23 l_v31 l_v32 l_v33 p_h p_v1 p_v2
# Fields to output
sampling.fields = velocity temperature
# Definitions for each probe
sampling/l_v11.type = LineSampler
sampling/l_v11.num_points = 128
sampling/l_v11.start = 66.666667 66.666667 3.125
sampling/l_v11.end = 66.666667 66.666667 396.875
sampling/l_v12.type = LineSampler
sampling/l_v12.num_points = 128
sampling/l_v12.start = 66.666667 200.0 3.125
sampling/l_v12.end = 66.666667 200.0 396.875
sampling/l_v13.type = LineSampler
sampling/l_v13.num_points = 128
sampling/l_v13.start = 66.666667 333.33333 3.125
sampling/l_v13.end = 66.666667 333.33333 396.875
sampling/l_v21.type = LineSampler
sampling/l_v21.num_points = 128
sampling/l_v21.start = 200.0 66.666667 3.125
sampling/l_v21.end = 200.0 66.666667 396.875
sampling/l_v22.type = LineSampler
sampling/l_v22.num_points = 128
sampling/l_v22.start = 200.0 200.0 3.125
sampling/l_v22.end = 200.0 200.0 396.875
sampling/l_v23.type = LineSampler
sampling/l_v23.num_points = 128
sampling/l_v23.start = 200.0 333.33333 3.125
sampling/l_v23.end = 200.0 333.33333 396.875
sampling/l_v31.type = LineSampler
sampling/l_v31.num_points = 128
sampling/l_v31.start = 333.33333 66.666667 3.125
sampling/l_v31.end = 333.33333 66.666667 396.875
sampling/l_v32.type = LineSampler
sampling/l_v32.num_points = 128
sampling/l_v32.start = 333.33333 200.0 3.125
sampling/l_v32.end = 333.33333 200.0 396.875
sampling/l_v33.type = LineSampler
sampling/l_v33.num_points = 128
sampling/l_v33.start = 333.33333 333.33333 3.125
sampling/l_v33.end = 333.33333 333.33333 396.875
sampling/p_h.type = PlaneSampler
sampling/p_h.axis1 = 1.0 0.0 0.0
sampling/p_h.axis2 = 0.0 1.0 0.0
sampling/p_h.origin = 0.0 0.0 15.0
sampling/p_h.num_points = 128 128
sampling/p_h.normal = 0.0 0.0 1.0
sampling/p_h.offsets = 0.0 15.0 35.0 85.0 185.0
sampling/p_v1.type = PlaneSampler
sampling/p_v1.axis1 = 0.0 0.0 1.0
sampling/p_v1.axis2 = 1.0 0.0 0.0
sampling/p_v1.origin = 0.0 200.0 0.0
sampling/p_v1.num_points = 128 128
sampling/p_v2.type = PlaneSampler
sampling/p_v2.axis1 = 0.0 0.0 1.0
sampling/p_v2.axis2 = 0.0 -1.0 0.0
sampling/p_v2.origin = 200.0 400.0 0.0
sampling/p_v2.num_points = 128 128

