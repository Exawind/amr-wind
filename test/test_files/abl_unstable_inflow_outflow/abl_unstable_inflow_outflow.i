#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =  4.0     # Max (simulated) time to evolve
time.max_step                =  8          # Max number of time steps
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.5        # Use this constant dt if > 0
time.cfl              =   0.95       # CFL factor
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
io.restart_file = "../abl_unstable_precursor/chk00002"
time.plot_interval            =  10       # Steps between plot files
time.checkpoint_interval      =  2       # Steps between checkpoint files
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity        =  0.0  0.0 -9.81  # Gravitational force (3D)
incflo.density        =  1.0             # Reference density
incflo.use_godunov = 1
incflo.godunov_type = "ppm_nolim"
transport.viscosity = 0.0
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = OneEqKsgsM84
incflo.physics = ABL
ICNS.source_terms = BoussinesqBuoyancy CoriolisForcing ABLForcing  # Uncomment for simulations without turbines
# ICNS.source_terms = BoussinesqBuoyancy CoriolisForcing ActuatorForcing BodyForce ABLMeanBoussinesq  # Uncomment for simulations with turbines
#--------- Additions by calc_inflow_stats.py ---------#
ABL.wall_shear_stress_type = "local"
ABL.inflow_outflow_mode = true
ABL.wf_velocity = 9.992795138134136 -0.0007138743053467138
ABL.wf_vmag = 10.016339397424394
ABL.wf_theta = 290.0178442735406
# BodyForce.magnitude = 0.0008901711909450327 0.0014550668111182673 0.0  # Uncomment for simulations with turbines
# BoussinesqBuoyancy.read_temperature_profile = true  # Uncomment for simulations with turbines
# BoussinesqBuoyancy.tprofile_filename = "../abl_unstable_precursor/avg_theta.dat"  # Uncomment for simulations with turbines
#-----------------------------------------------------#
TKE.source_terms = KsgsM84Src
# TKE.interpolation="PiecewiseConstant"  # Use if amr.max_level > 1
BoussinesqBuoyancy.reference_temperature = 290.0
CoriolisForcing.east_vector = 1.0 0.0 0.0
CoriolisForcing.north_vector = 0.0 1.0 0.0
CoriolisForcing.latitude = 90.0
ABLForcing.abl_forcing_height= 100.0
incflo.velocity = 10.0 0.0 0.0
ABL.reference_temperature = 290.0
ABL.temperature_heights = 0.0 1000.0 1100.0 2100.0
ABL.temperature_values = 290.0 290.0 298.0 301.0
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
ABL.stats_output_frequency = 1
ABL.bndry_file = "../abl_unstable_precursor/bndry_files"
ABL.bndry_io_mode = 1
ABL.bndry_planes = xlo
ABL.bndry_output_start_time = 1.0
ABL.bndry_var_names = velocity temperature tke

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 36 36 24    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy
amr.blocking_factor     = 4
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.       0.     0.  # Lo corner coordinates
geometry.prob_hi        =   3000.  3000.  2000.  # Hi corner coordinates
geometry.is_periodic    =   0   1   0   # Periodicity x y z (0/1)
# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.225
xlo.temperature = 290.0
xlo.tke = 0
xhi.type = "pressure_outflow"
zlo.type =   "wall_model"
zlo.temperature_type = "wall_model"
zhi.type =   "slip_wall"
zhi.temperature_type = "fixed_gradient"
zhi.temperature = 0.003
zlo.tke_type = "fixed_gradient"
incflo.verbose          =   0
