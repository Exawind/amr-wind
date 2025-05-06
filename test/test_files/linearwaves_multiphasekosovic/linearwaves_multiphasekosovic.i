#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION CONTROL         #
#.......................................#
time.stop_time                           = 7200
time.max_step                            = -1          # Max number of time steps
time.fixed_dt                            = -1          # Use this constant dt if > 0
time.initial_dt                          = 0.05
time.cfl                                 = 0.9         # CFL factor

time.plot_interval                       = 1000       # Steps between plot files
time.checkpoint_interval                 = 25000       # Steps between checkpoint files
io.output_default_variables              = false
io.outputs                               = velocity terrain_blank  ow_velocity
incflo.use_godunov                       = 1
incflo.diffusion_type                    = 1
incflo.godunov_type                      = weno_z
  incflo.mflux_type                        = upwind
turbulence.model                         = MultiPhaseKosovic
incflo.gravity                           = 0.  0. -9.81  # Gravitational force (3D)
  transport.model                          = TwoPhaseTransport
transport.viscosity_fluid1               = 1e-3
transport.viscosity_fluid2               = 1e-5
transport.laminar_prandtl_fluid1         = 7.2
transport.laminar_prandtl_fluid2         = 0.7
transport.turbulent_prandtl              = 0.3333
MultiPhase.density_fluid1                = 1020.
MultiPhase.density_fluid2                = 1.225
MultiPhase.water_level                   = 0.
incflo.verbose                           =   2          # incflo_level

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            GEOMETRY & BCs             #
#.......................................#
geometry.prob_lo                         = 0.       0.    -200.  # Lo corner coordinates
geometry.prob_hi                         = 3000 3000 1600  # Hi corner coordinates; 600 / 128 = 4.6875 m resolution
amr.n_cell                               = 384 384 904
amr.max_level                            = 0
geometry.is_periodic                     = 1   1   0  
zlo.type                                 = slip_wall
zlo.temperature_type                     = fixed_gradient
zlo.temperature                          = 0.0										 
zhi.type                                 = slip_wall
zhi.temperature_type                     = fixed_gradient
zhi.temperature                          = 0.003

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.physics                           = MultiPhase ABL OceanWaves
ICNS.source_terms                        = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing GravityForcing NonLinearSGSTerm
ICNS.use_perturb_pressure                = true
incflo.velocity                          = 10.0 0.0 0.0
GeostrophicForcing.geostrophic_wind      = 10.0 0.0 0.0
GeostrophicForcing.wind_forcing_off_height  = 12.0
GeostrophicForcing.wind_forcing_ramp_height = 16.0
CoriolisForcing.latitude                 = 90.0
CoriolisForcing.north_vector             = 0.0 1.0 0.0
CoriolisForcing.east_vector              = 1.0 0.0 0.0
BoussinesqBuoyancy.reference_temperature = 290
RayleighDamping.reference_velocity       = 10.0 0 0
RayleighDamping.length_sloped_damping    = 400
RayleighDamping.length_complete_damping  = 200
RayleighDamping.time_scale               = 20.0
RayleighDamping.force_coord_directions   = 0 0 1								
ABL.reference_temperature                = 290
ABL.temperature_heights                  = 0.0 300 1300 2300
ABL.temperature_values                   = 290 290 293 296
ABL.perturb_temperature                  = false
ABL.cutoff_height                        = 50.0
ABL.perturb_velocity                     = true
ABL.perturb_ref_height                   = 50.0
ABL.Uperiods                             = 4.0
ABL.Vperiods                             = 4.0
ABL.deltaU                               = 1e-5
ABL.deltaV                               = 1e-5
ABL.kappa                                = .41
ABL.surface_roughness_z0                 = 0.0001
ABL.surface_temp_flux                    = 0.0
ABL.wall_shear_stress_type               = local
											       
OceanWaves.label                         = Wave1
OceanWaves.Wave1.type                    = LinearWaves
OceanWaves.Wave1.wave_height             = 6
OceanWaves.Wave1.wave_length             = 400.0
OceanWaves.Wave1.water_depth             = 200.0
OceanWaves.Wave1.relax_zone_gen_length   = 400.0
OceanWaves.Wave1.relax_zone_out_length   = 800.0
OceanWaves.Wave1.initialize_wave_field   = true

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          POST-Processing              #
#.......................................#


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              AVERAGING                #
#.......................................#


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            MESH REFINEMENT            #
#.......................................#
tagging.labels                           = static
tagging.static.type                      = CartBoxRefinement
tagging.static.static_refinement_def     = static_box.txt

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               TURBINES                #
#.......................................#
mac_proj.num_pre_smooth                            = 8
mac_proj.num_post_smooth                           = 8
mac_proj.mg_rtol                                   = -1
mac_proj.mg_atol                                   = 1e-4
mac_proj.maxiter                                   = 25
mac_proj.fmg_maxiter                               = 4
nodal_proj.num_pre_smooth                          = 8
nodal_proj.num_post_smooth                         = 8
nodal_proj.mg_rtol                                 = -1
nodal_proj.mg_atol                                 = 1e-4
nodal_proj.maxiter                                 = 25
nodal_proj.fmg_maxiter                             = 4
diffusion.mg_rtol                                  = -1
diffusion.mg_atol                                  = 1e-4
temperature_diffusion.mg_rtol                      = -1
temperature_diffusion.mg_atol                      = 1e-4

