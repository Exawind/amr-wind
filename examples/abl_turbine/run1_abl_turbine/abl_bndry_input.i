#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   20099.9     # Max (simulated) time to evolve
time.max_step                =   60000          # Max number of time steps
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.1        # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
io.restart_file = "../amrwind_abl/chk40000"
#io.restart_file = "chk40002"
time.plot_interval            =  500       # Steps between plot files
time.checkpoint_interval      =  500       # Steps between checkpoint files

#time.regrid_interval = 1

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity          =   0.  0. -9.81  # Gravitational force (3D)
incflo.density          = 1.0          # Reference density
incflo.use_godunov = 1
transport.viscosity = 1.0e-5
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
#turbulence.model = Smagorinsky
turbulence.model = OneEqKsgsM84
#TKE.source_terms = KsgsM84Src
Smagorinsky_coeffs.Cs = 0.135
incflo.physics = ABL Actuator
#incflo.physics = ABL
ICNS.source_terms = CoriolisForcing GeostrophicForcing ActuatorForcing
#ICNS.source_terms = CoriolisForcing GeostrophicForcing
BoussinesqBuoyancy.reference_temperature = 290.0
ABL.reference_temperature = 290.0
CoriolisForcing.east_vector = 1.0 0.0 0.0
CoriolisForcing.north_vector = 0.0 1.0 0.0
CoriolisForcing.latitude = 90.0
CoriolisForcing.rotational_time_period = 125663.706143592
GeostrophicForcing.geostrophic_wind = 10.0 0.0 0.0
incflo.velocity = 10.0 0.0 0.0
ABL.temperature_heights = 0.0 2000.0
ABL.temperature_values = 290.0 290.0
ABL.perturb_temperature = false
ABL.cutoff_height = 50.0
ABL.perturb_velocity = true
ABL.perturb_ref_height = 50.0
ABL.Uperiods = 4.0
ABL.Vperiods = 4.0
ABL.deltaU = 1.0
ABL.deltaV = 1.0
ABL.kappa = .41
ABL.surface_roughness_z0 = 0.01
ABL.bndry_file = "../amrwind_abl/bndry_file.nc"
ABL.bndry_io_mode = 1
ABL.bndry_var_names = velocity temperature tke
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 512 512 128    # Grid cells at coarsest AMRlevel
#amr.n_cell              = 128 128 32    # Grid cells at coarsest AMRlevel
amr.max_level           = 2           # Max AMR level in hierarchy 
#amr.max_level           = 1           # Max AMR level in hierarchy 
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.       0.     0.  # Lo corner coordinates
geometry.prob_hi        =   5000.  5000.  1000.  # Hi corner coordinates
geometry.is_periodic    =   0   1   0   # Periodicity x y z (0/1)
incflo.delp             =   0.  0.  0.  # Prescribed (cyclic) pressure gradient


# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.0
xlo.temperature = 0.0
xhi.type = "pressure_outflow"
xlo.tke = 0
#ylo.type = "mass_inflow"
#ylo.density = 1.0
#ylo.temperature = 0.0
#yhi.type = "pressure_outflow"
zlo.type =   "wall_model"
zhi.type =   "slip_wall"
zhi.temperature_type = "fixed_gradient"
zhi.temperature = 0.0
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              VERBOSITY                #
#.......................................#
incflo.verbose          =   0          # incflo_level


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        Mesh refinement                #
#.......................................#
tagging.labels = refine0
tagging.refine0.type = GeometryRefinement
tagging.refine0.shapes = c0 c1

tagging.refine0.c0.type = cylinder
tagging.refine0.c0.start = 1450.0 2500. 90.
tagging.refine0.c0.end = 1600.0 2500. 90.
tagging.refine0.c0.outer_radius = 189.0

tagging.refine0.c1.type = cylinder
tagging.refine0.c1.start = 1500.0 2500. 90.
#tagging.refine0.c1.end = 2000.0 2500. 90.
tagging.refine0.c1.end = 1550.0 2500. 90.
tagging.refine0.c1.outer_radius = 157.5


Actuator.labels = WTG01
Actuator.type = TurbineFastLine

Actuator.TurbineFastLine.rotor_diameter = 126.0
Actuator.TurbineFastLine.hub_height = 90.0
Actuator.TurbineFastLine.num_points_blade = 64
Actuator.TurbineFastLine.num_points_tower = 12
Actuator.TurbineFastLine.epsilon = 5.0 5.0 5.0
Actuator.TurbineFastLine.epsilon_tower = 5.0 5.0 5.0
Actuator.TurbineFastLine.openfast_start_time = 0.0
Actuator.TurbineFastLine.openfast_stop_time = 2000.0
Actuator.TurbineFastLine.nacelle_drag_coeff = 0.0
Actuator.TurbineFastLine.nacelle_area = 0.0
Actuator.TurbineFastLine.output_frequency = 10
Actuator.TurbineFastLine.density = 1.225

Actuator.WTG01.base_position = 1500.0 2500. 0.0
Actuator.WTG01.openfast_input_file = "fast_inp/nrel5mw.fst"


