#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   22000.0     # Max (simulated) time to evolve
time.max_step                =   10          # Max number of time steps
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.4        # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
io.KE_int = 1
io.restart_file = "../abl_bndry_output/chk00005"
time.plot_interval            =  1       # Steps between plot files
time.checkpoint_interval      =  1       # Steps between checkpoint files
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity          =   0.  0. -9.81  # Gravitational force (3D)
incflo.density          = 1.0          # Reference density
incflo.use_godunov = 1
transport.viscosity = 1.0e-5
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Smagorinsky
Smagorinsky_coeffs.Cs = 0.135
incflo.physics = ABL
ICNS.source_terms = CoriolisForcing GeostrophicForcing
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
ABL.bndry_file = "../abl_bndry_output/bndry_file.nc"
ABL.bndry_io_mode = 1
ABL.bndry_var_names = velocity temperature
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 48 48 48    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy 
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.       0.     0.  # Lo corner coordinates
geometry.prob_hi        =   1000.  1000.  1000.  # Hi corner coordinates
geometry.is_periodic    =   0   0   0   # Periodicity x y z (0/1)
incflo.delp             =   0.  0.  0.  # Prescribed (cyclic) pressure gradient
# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.0
xlo.temperature = 0.0
xhi.type = "pressure_outflow"
ylo.type = "mass_inflow"
ylo.density = 1.0
ylo.temperature = 0.0
yhi.type = "pressure_outflow"
zlo.type =   "wall_model"
zhi.type =   "slip_wall"
zhi.temperature_type = "fixed_gradient"
zhi.temperature = 0.0
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              VERBOSITY                #
#.......................................#
incflo.verbose          =   0          # incflo_level
