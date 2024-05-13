.. _turbine:

Turbine Simulation Walkthrough
==============================

Now that we have run our precursor simulation and saved inflow boundary condition data, we can run a turbine simulation. Here is the config file:
```
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION CONTROL         #
#.......................................#
time.stop_time                           = 7801.0             # Max (simulated) time to evolve [s]
time.max_step                            = -1          # Max number of time steps; -1 means termination set by timestamps
time.fixed_dt                            = 0.125        # Use this constant dt if > 0
time.cfl                                 = 0.95         # CFL factor

time.plot_interval                       = 1200       # Steps between plot files
time.checkpoint_interval                 = 1200       # Steps between checkpoint files
ABL.bndry_file                           = ../02_atmosphere/precursor/bndry_file.native
ABL.bndry_io_mode                        = 1          # 0 = write, 1 = read
ABL.bndry_planes                         = ylo xlo
ABL.bndry_output_start_time              = 7200.0
ABL.bndry_var_names                      = velocity temperature tke

incflo.physics                           = ABL Actuator
io.restart_file                          = ../02_atmosphere/spinup/chk14400   
incflo.use_godunov                       = 1
incflo.godunov_type                      = weno_z                 
turbulence.model                         = OneEqKsgsM84  # For neutral ABL, use "Smagorinsky"
TKE.source_terms                         = KsgsM84Src
TKE.interpolation                        = PiecewiseConstant          
incflo.gravity                           = 0.  0. -9.81  # Gravitational force (3D)
incflo.density                           = 1.225          # Reference density; make sure this agrees with OpenFAST values
transport.viscosity                      = 1.0e-5
transport.laminar_prandtl                = 0.7
transport.turbulent_prandtl              = 0.3333

incflo.verbose                           =   0          # incflo_level

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            GEOMETRY & BCs             #
#.......................................#
geometry.prob_lo                         = 0.       0.     0.  # Lo corner coordinates
geometry.prob_hi                         = 2560.  2560.  1280.  # Hi corner coordinates
amr.n_cell                               = 128 128 64    # Grid cells at coarsest AMRlevel
amr.max_level                            = 2           # Max AMR level in hierarchy 
geometry.is_periodic                     = 0   0   0   # Periodicity x y z (0/1)
incflo.delp                              = 0.  0.  0.  # Prescribed (cyclic) pressure gradient

xlo.type                                 = mass_inflow         
xlo.density                              = 1.225               
xlo.temperature                          = 290.0               
xlo.tke                                  = 0.0
xhi.type                                 = pressure_outflow    

ylo.type                                 = mass_inflow         
ylo.density                              = 1.225               
ylo.temperature                          = 290.0               
ylo.tke                                  = 0.0
yhi.type                                 = pressure_outflow     

zlo.type                                 = wall_model
zhi.type                                 = slip_wall
zhi.temperature_type                     = fixed_gradient
zhi.temperature                          = 0.003

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
ICNS.source_terms                        = BoussinesqBuoyancy CoriolisForcing BodyForce ABLMeanBoussinesq ActuatorForcing
##--------- Additions by calc_inflow_stats.py ---------#
ABL.wall_shear_stress_type = "local"
ABL.inflow_outflow_mode = true
ABL.wf_velocity = 7.612874161922772 0.05167978610261268
ABL.wf_vmag = 7.655337919472146
ABL.wf_theta = 291.0187046202964
BodyForce.magnitude = 0.00034839850789284793 0.0009712494385077595 0.0
BoussinesqBuoyancy.read_temperature_profile = true
BoussinesqBuoyancy.tprofile_filename = avg_theta.dat
##-----------------------------------------------------#
incflo.velocity                          = 10.0 0.0 0.0
ABLForcing.abl_forcing_height            = 86.5
CoriolisForcing.latitude                 = 36.607322      # Southern Great Planes
CoriolisForcing.north_vector             = 0.0 1.0 0.0
CoriolisForcing.east_vector              = 1.0 0.0 0.0
BoussinesqBuoyancy.reference_temperature = 290.0
ABL.reference_temperature                = 290.0
ABL.temperature_heights                  = 0.0 600.0 700.0 1700.0    # Make sure the top height here goes above the domain height
ABL.temperature_values                   = 290.0 290.0 298.0 301.0
ABL.perturb_temperature                  = true
ABL.cutoff_height                        = 50.0
ABL.perturb_velocity                     = true
ABL.perturb_ref_height                   = 50.0
ABL.Uperiods                             = 4.0
ABL.Vperiods                             = 4.0
ABL.deltaU                               = 1.0
ABL.deltaV                               = 1.0
ABL.kappa                                = .40
ABL.surface_roughness_z0                 = 0.01  # meters
ABL.surface_temp_flux                    = 0.05  # Surface temperature flux [K-m/s]

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          POST-Processing              #
#.......................................#
#io.output_hdf5_plotfile                  = true  # Uncomment these two lines if save to .hdf instead of plt#####/
#io.hdf5_compression                      = "ZFP_ACCURACY@0.001"

incflo.post_processing                   = sampling averaging

# --- Sampling parameters ---
sampling.output_frequency                = 8
sampling.fields                          = velocity temperature

#---- sample defs ----
sampling.labels                          = xy-domain xz-domain 

sampling.xy-domain.type                  = PlaneSampler        
sampling.xy-domain.num_points            = 256 256             
sampling.xy-domain.origin                = 0.0 0.0 91.0      
sampling.xy-domain.axis1                 = 2550.0 0.0 0.0      
sampling.xy-domain.axis2                 = 0.0 2550.0 0.0      
sampling.xy-domain.normal                = 0.0 0.0 1.0         
sampling.xy-domain.offsets               = -63.45 0.0 63.45  

sampling.xz-domain.type                  = PlaneSampler        
sampling.xz-domain.num_points            = 256 128              
sampling.xz-domain.origin                = 0.0 1280.0 0.0         
sampling.xz-domain.axis1                 = 2550.0 0.0 0.0      
sampling.xz-domain.axis2                 = 0.0 0.0 1270.0   

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              AVERAGING                #
#.......................................#
averaging.type                           = TimeAveraging
averaging.labels                         = means stress

averaging.averaging_window               = 60.0
averaging.averaging_start_time           = 7200.0

averaging.means.fields                   = velocity
averaging.means.averaging_type           = ReAveraging

averaging.stress.fields                  = velocity
averaging.stress.averaging_type          = ReynoldsStress

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            MESH REFINEMENT            #
#.......................................#
tagging.labels                           = T0_level_0_zone T1_level_0_zone T2_level_0_zone T0_level_1_zone T1_level_1_zone T2_level_1_zone

# 1st refinement level
tagging.T0_level_0_zone.type             = GeometryRefinement  
tagging.T0_level_0_zone.shapes           = T0_level_0_zone     
tagging.T0_level_0_zone.level            = 0                   
tagging.T0_level_0_zone.T0_level_0_zone.type = box                 
tagging.T0_level_0_zone.T0_level_0_zone.origin = 520.0 1040.0 0.0  # -1D, -2D
tagging.T0_level_0_zone.T0_level_0_zone.xaxis = 360.0 0.0 0.0
tagging.T0_level_0_zone.T0_level_0_zone.yaxis = 0.0 480.0 0.0
tagging.T0_level_0_zone.T0_level_0_zone.zaxis = 0.0 0.0 360.0

tagging.T1_level_0_zone.type             = GeometryRefinement  
tagging.T1_level_0_zone.shapes           = T1_level_0_zone     
tagging.T1_level_0_zone.level            = 0                   
tagging.T1_level_0_zone.T1_level_0_zone.type = box                 
tagging.T1_level_0_zone.T1_level_0_zone.origin = 1160.0 1040.0 0.0  # -1D, -2D
tagging.T1_level_0_zone.T1_level_0_zone.xaxis = 360.0 0.0 0.0
tagging.T1_level_0_zone.T1_level_0_zone.yaxis = 0.0 480.0 0.0
tagging.T1_level_0_zone.T1_level_0_zone.zaxis = 0.0 0.0 360.0

tagging.T2_level_0_zone.type             = GeometryRefinement  
tagging.T2_level_0_zone.shapes           = T2_level_0_zone     
tagging.T2_level_0_zone.level            = 0                   
tagging.T2_level_0_zone.T2_level_0_zone.type = box                 
tagging.T2_level_0_zone.T2_level_0_zone.origin = 1800.0 1040.0 0.0  # -1D, -2D
tagging.T2_level_0_zone.T2_level_0_zone.xaxis = 360.0 0.0 0.0
tagging.T2_level_0_zone.T2_level_0_zone.yaxis = 0.0 480.0 0.0
tagging.T2_level_0_zone.T2_level_0_zone.zaxis = 0.0 0.0 360.0

# 2nd refinement level
tagging.T0_level_1_zone.type             = GeometryRefinement  
tagging.T0_level_1_zone.shapes           = T0_level_1_zone     
tagging.T0_level_1_zone.level            = 1                   
tagging.T0_level_1_zone.T0_level_1_zone.type = box                 
tagging.T0_level_1_zone.T0_level_1_zone.origin = 580.0 1100.0 20.0  # -0.5D, -1.5D
tagging.T0_level_1_zone.T0_level_1_zone.xaxis = 180.0 0.0 0.0
tagging.T0_level_1_zone.T0_level_1_zone.yaxis = 0.0 360.0 0.0
tagging.T0_level_1_zone.T0_level_1_zone.zaxis = 0.0 0.0 180.0

tagging.T1_level_1_zone.type             = GeometryRefinement  
tagging.T1_level_1_zone.shapes           = T1_level_1_zone     
tagging.T1_level_1_zone.level            = 1                   
tagging.T1_level_1_zone.T1_level_1_zone.type = box                 
tagging.T1_level_1_zone.T1_level_1_zone.origin = 1220.0 1100.0 20.0  # -0.5D, -1.5D
tagging.T1_level_1_zone.T1_level_1_zone.xaxis = 180.0 0.0 0.0
tagging.T1_level_1_zone.T1_level_1_zone.yaxis = 0.0 360.0 0.0
tagging.T1_level_1_zone.T1_level_1_zone.zaxis = 0.0 0.0 180.0

tagging.T2_level_1_zone.type             = GeometryRefinement  
tagging.T2_level_1_zone.shapes           = T2_level_1_zone     
tagging.T2_level_1_zone.level            = 1                   
tagging.T2_level_1_zone.T2_level_1_zone.type = box                 
tagging.T2_level_1_zone.T2_level_1_zone.origin = 1860.0 1100.0 20.0  # -0.5D, -1.5D
tagging.T2_level_1_zone.T2_level_1_zone.xaxis = 180.0 0.0 0.0
tagging.T2_level_1_zone.T2_level_1_zone.yaxis = 0.0 360.0 0.0
tagging.T2_level_1_zone.T2_level_1_zone.zaxis = 0.0 0.0 180.0

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               TURBINES                #
#.......................................#
Actuator.labels                          = T0 T1 T2       

Actuator.T0.type                         = TurbineFastDisk     
Actuator.T0.openfast_input_file          = T0_AMRWind_AWAKEN/NREL-2p8-127.fst
Actuator.T0.base_position                = 640.0 1280.0 0.0   
Actuator.T0.rotor_diameter               = 126.9               
Actuator.T0.hub_height                   = 86.5                
Actuator.T0.num_points_blade             = 64                  
Actuator.T0.num_points_tower             = 12                  
Actuator.T0.epsilon                      = 5.0 5.0 5.0         
Actuator.T0.epsilon_tower                = 5.0 5.0 5.0         
Actuator.T0.openfast_start_time          = 0.0                 
Actuator.T0.openfast_stop_time           = 99999.0              
Actuator.T0.nacelle_drag_coeff           = 0.0                 
Actuator.T0.nacelle_area                 = 0.0                 
Actuator.T0.yaw                          = 0.0   
Actuator.T0.output_frequency             = 10            

Actuator.T1.type                         = TurbineFastDisk     
Actuator.T1.openfast_input_file          = T1_AMRWind_AWAKEN/NREL-2p8-127.fst
Actuator.T1.base_position                = 1280.0 1280.0 0.0   
Actuator.T1.rotor_diameter               = 126.9               
Actuator.T1.hub_height                   = 86.5                
Actuator.T1.num_points_blade             = 64                  
Actuator.T1.num_points_tower             = 12                  
Actuator.T1.epsilon                      = 5.0 5.0 5.0         
Actuator.T1.epsilon_tower                = 5.0 5.0 5.0         
Actuator.T1.openfast_start_time          = 0.0                 
Actuator.T1.openfast_stop_time           = 99999.0                
Actuator.T1.nacelle_drag_coeff           = 0.0                 
Actuator.T1.nacelle_area                 = 0.0                 
Actuator.T1.yaw                          = 0.0   
Actuator.T1.output_frequency             = 10              

Actuator.T2.type                         = TurbineFastDisk     
Actuator.T2.openfast_input_file          = T2_AMRWind_AWAKEN/NREL-2p8-127.fst
Actuator.T2.base_position                = 1920.0 1280.0 0.0   
Actuator.T2.rotor_diameter               = 126.9               
Actuator.T2.hub_height                   = 86.5                
Actuator.T2.num_points_blade             = 64                  
Actuator.T2.num_points_tower             = 12                  
Actuator.T2.epsilon                      = 5.0 5.0 5.0         
Actuator.T2.epsilon_tower                = 5.0 5.0 5.0         
Actuator.T2.openfast_start_time          = 0.0                 
Actuator.T2.openfast_stop_time           = 99999.0                
Actuator.T2.nacelle_drag_coeff           = 0.0                 
Actuator.T2.nacelle_area                 = 0.0                 
Actuator.T2.yaw                          = 0.0   
Actuator.T2.output_frequency             = 10         

```

This file looks like the precursor config file, except the following changes:
* We're now reading in boundary condition data, not writing it out. Similarly, the x- and y- boundaries are no longer periodic, and we specify some extra characteristics about the inflow.
* `incflo.physics` now includes `Actuator`
* Similarly, we drop `ABLForcing` from `ICNS.source_terms`, and we replace it with three new forces: `BodyForce ABLMeanBoussinesq ActuatorForcing`. We then provide extra information about these forces. More on how to calculate these values down below.
* We add information about the turbines into the config file

When running an inflow-outflow simulation, you need to calculate information about `BodyForce` and `ABLMeanBoussinesq` from the precursor simulation. Do do this for this example, run the `calc_inflowoutflow_stats.py` as follow:
```
python calc_inflowoutflow_stats.py -sf /scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/precursor/post_processing/abl_statistics14400.nc -ts 7200 -te 9000
```

This script gives you two important things:
* Information that you should copy-paste into the config file
* An `avg_theta.dat` file with important information for `ABLMeanBoussinesq`. It is critical that you copy or symbolically link this file into the same place as where you config file will run (usually /scratch/), otherwise AMR-Wind will fail with a mysterious error message

Because we're simulating turbines, we also need to include OpenFAST files for each of the turbines. In this demo, we use 2.8 MW turbines from NREL's open source [turbine repo](https://github.com/NREL/openfast-turbine-models/tree/master/IEA-scaled/NREL-2.8-127/OpenFAST). When simulating OpenFAST turbines through AMR-Wind instead of directly through OpenFAST, it is important to make the follow changes to the OpenFAST files:
* AeroDyn: Make sure `WakeMod` is 0
* ElastoDyn: Set the initial RPM `RotSpeed` and inital yaw angle `NacYaw` to reasonable values
* `*.fst`: Set `CompInflow` to be 2 and `OutFileFmt` to be 1
* ServoDyn: Make sure `DLL_FileName` points to a `libdiscon.so` file from ROSCO

When all is said and done, the /scratch/ directory where I run the turbine simulation looks like
```
amr_wind  avg_theta.dat	T0_AMRWind_AWAKEN/  T1_AMRWind_AWAKEN/	T2_AMRWind_AWAKEN/  turbine.i