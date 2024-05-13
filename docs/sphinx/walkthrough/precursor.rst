.. _precursor:

Precursor (ABL) Walkthrough
=====================

Now that we have AMR-Wind compiled, we can start to run simulations. In this section, we'll run two simulations of a weakly convective atmospheric boundary layer. The first simulation will be of a transient "spinup" phase, and the second simulation will serve as a "precursor" for a turbine simulation.

Before we run anything, we first need to set up the [input](https://exawind.github.io/amr-wind/user/inputs.html) file (aka the "configuration"/"config" file). This is a text file, and its filename traditionally ends with `.i` or `.inp`. There are two general approaches to set one up: manually, or with the use of [amrwind-frontend](https://github.com/lawrenceccheung/amrwind-frontend). If I am setting up a big simulation, with many turbines and many refinement zones, I use amrwind-frontend. It shows the locations of all those objects, which is a great sanity check for expensive simulations. If I am running a simpler simulation, I will copy-paste text from a simulation that I have successfully run in the past. 

### Spinup
Here is the content of our spinup config file:
```
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION CONTROL         #
#.......................................#
time.stop_time                           = 7201.0             # Max (simulated) time to evolve [s]
time.max_step                            = -1          # Max number of time steps; -1 means termination set by timestamps
time.fixed_dt                            = 0.5        # Use this constant dt if > 0
time.cfl                                 = 0.95         # CFL factor

time.plot_interval                       = 1800       # Steps between plot files
time.checkpoint_interval                 = 1800       # Steps between checkpoint files
#ABL.bndry_file                           = bndry_file.native
#ABL.bndry_io_mode                        = 0          # 0 = write, 1 = read
#ABL.bndry_planes                         = ylo xlo
#ABL.bndry_output_start_time              = 7200.0
#ABL.bndry_var_names                      = velocity temperature tke

incflo.physics                           = ABL # Actuator
#io.restart_file                          = ../spinup/chk07200   
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
amr.max_level                            = 0           # Max AMR level in hierarchy 
geometry.is_periodic                     = 1   1   0   # Periodicity x y z (0/1)
incflo.delp                              = 0.  0.  0.  # Prescribed (cyclic) pressure gradient

#xlo.type                                 = mass_inflow         
#xlo.density                              = 1.225               
#xlo.temperature                          = 290.0               
#xlo.tke                                  = 0.0
#xhi.type                                 = pressure_outflow    

#ylo.type                                 = mass_inflow         
#ylo.density                              = 1.225               
#ylo.temperature                          = 290.0               
#ylo.tke                                  = 0.0
#yhi.type                                 = pressure_outflow     

zlo.type                                 = wall_model
zhi.type                                 = slip_wall
zhi.temperature_type                     = fixed_gradient
zhi.temperature                          = 0.003

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
ICNS.source_terms                        = BoussinesqBuoyancy CoriolisForcing ABLForcing
##--------- Additions by calc_inflow_stats.py ---------#
#ABL.wall_shear_stress_type = "local"
#ABL.inflow_outflow_mode = true
#ABL.wf_velocity = XX XX
#ABL.wf_vmag = XX
#ABL.wf_theta = XX
#BodyForce.magnitude = XX XX XX
#BoussinesqBuoyancy.read_temperature_profile = true
#BoussinesqBuoyancy.tprofile_filename = avg_theta.dat
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

incflo.post_processing                   = sampling # averaging

# --- Sampling parameters ---
sampling.output_frequency                = 900                 
sampling.fields                          = velocity temperature

#---- sample defs ----
sampling.labels                          = xy-domain xz-domain 

sampling.xy-domain.type                  = PlaneSampler        
sampling.xy-domain.num_points            = 256 256             
sampling.xy-domain.origin                = 0.0 0.0 86.5      
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
#averaging.type                           = TimeAveraging
#averaging.labels                         = means stress

#averaging.averaging_window               = 600.0
#averaging.averaging_start_time           = 7200.0

#averaging.means.fields                   = velocity
#averaging.means.averaging_type           = ReAveraging

#averaging.stress.fields                  = velocity
#averaging.stress.averaging_type          = ReynoldsStress

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            MESH REFINEMENT            #
#.......................................#


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               TURBINES                #
#.......................................#
```
To get an idea of what each of these lines means, see [here](https://exawind.github.io/amr-wind/user/inputs.html). Sometimes new input variables are added but aren't located in the docs, and in that case, you will need to look into the guts of the AMR-Wind code to deduce what those variables do.

We are going to run the spinup simulation for a duration of two hours. In this tutorial, we are simulation a fairly small domain with a coarse uniform grid resolution of 20 m. The spinup simulation does not employ any grid refinement.

I am running this simulation on NREL's supercomputer Eagle, and I have this file placed in the directory `/projects/wakedynamics/orybchuk/amr-wind-tutorial/02_atmosphere/spinup`. The contents of this directory are 
```
analyze_spinup.ipynb  logs/  runfile.sbatch  spinup.i
``` 

I like executing simulations in Eagle's scratch directory, so I also created `/scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/spinup`. I `cd` into that `/scratch/` directory and then symbolically link to the config file in projects, with the command
```
ln -sf /projects/wakedynamics/orybchuk/amr-wind-tutorial/02_atmosphere/spinup/spinup.i .
```

Now, the contents of the scratch directory are 
```
spinup.i
```

To kick off the job on the Eagle's compute nodes, I use `runfile.sbatch`. Its contents are
```
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --account=awaken
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XX@gmail.com
#SBATCH --output=logs/job_output_filename.%j.out  # %j will be replaced with the job ID
#SBATCH --partition=debug
#SBATCH --time=00:59:00
# #SBATCH --partition=short
# #SBATCH --time=03:59:00
# #SBATCH --partition=standard
# #SBATCH --time=47:59:00

module purge
module load gcc/8.4.0 mpt mkl cmake

export EXAWIND_DIR=/nopt/nrel/ecom/exawind/exawind-2020-09-21/install/gcc
export MPI_TYPE_DEPTH=15
export MPI_IB_CONGESTED=true
export MPI_XPMEM_ENABLED=disabled


cd /scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/spinup

rm -rf post_processing
ln -sf /projects/awaken/orybchuk/spack-june22/amr-wind/spack-build-4ixvlaf/amr_wind .
srun -n 72 -c 1 --cpu_bind=cores amr_wind spinup.i
```

In this script, I symbolically link my AMR-Wind executable into the /scratch/ directory. I also run on two nodes, each of which have 36 cores, hence `-n 72`. Some of the other variables in this script, like `EXAWIND_DIR`, are specific to Eagle. You might not need to specify these in order to get AMR-Wind to run.

To kick off the job, go into the /projects/ directory, and execute
```
sbatch runfile.sbatch
```

You can view the job log by opening up `logs/job_output_filename.###.out`. Note: you need to manually create `logs/` before kicking off your simulation, otherwise your job will silently fail because of SLURM issues.

If your job is successful, the log file will look something like this at the end:
```
Writing plot file       plt14402 at time 7201
Writing checkpoint file chk14402 at time 7201
Time spent in InitData():    0.816465397
Time spent in Evolve():      2404.252191
Unused ParmParse Variables:
  [TOP]::incflo.delp(nvals = 3)  :: [0., 0., 0.]

AMReX (22.05-29-g1305eb3d364d) finalized
```

And the /scratch/ directory will now look something like
```
amr_wind@   chk05400/	chk09000/  chk12600/  chk14402/   core.27800  core.9246   plt05400/  plt09000/  plt12600/  plt14402/	      spinup.i@
chk00000/  chk01800/		    chk03600/		     chk07200/	chk10800/  chk14400/  core.27159  core.8529   plt00000/  plt01800/		  plt03600/		   plt07200/  plt10800/  plt14400/  post_processing/
```

Once the spinup simulation is done, it is important to sanity check that the fields make sense. One quick test to do this would be to look at the evolution of horizontally averaged vertical profiles. The Jupyter notebook `analyze_spinup.ipynb` contains code to help with this.

You can also open the `plt#####` files using Paraview or other software that AMReX is compatible with. These files show a volume of the instantaneous fields at that timestep.

### Precursor simulation
After sufficiently spinning up turbulence, I kick off a "precursor simulation". 

In the context of wind turbine LES, a precursor is a simulation that is run without a turbine for the explicit purpose of generating inflow boundary conditions. Due to the way LES works, spinup and precursor simulations are almost always run with cyclic boundary conditions. This means wind that exits the outflow simulation is then recirculated back into the inflow. If we want to simulate a statistically homogeneous atmosphere, that's fine. However, this is problematic if you have a wind turbine---turbines generate wakes, and we don't want wakes recirculating back into the inflow. So we run wind turbine simulations with a prescribed "inflow boundary condition" (where the wind data comes from the precursor) and an "outflow boundary condition" (usually a pressure BC).

Here is the content of the precursor simulation.
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
ABL.bndry_file                           = bndry_file.native
ABL.bndry_io_mode                        = 0          # 0 = write, 1 = read
ABL.bndry_planes                         = ylo xlo
ABL.bndry_output_start_time              = 7200.0
ABL.bndry_var_names                      = velocity temperature tke

incflo.physics                           = ABL # Actuator
io.restart_file                          = ../spinup/chk14400   
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
geometry.is_periodic                     = 1   1   0   # Periodicity x y z (0/1)
incflo.delp                              = 0.  0.  0.  # Prescribed (cyclic) pressure gradient

#xlo.type                                 = mass_inflow         
#xlo.density                              = 1.225               
#xlo.temperature                          = 290.0               
#xlo.tke                                  = 0.0
#xhi.type                                 = pressure_outflow    

#ylo.type                                 = mass_inflow         
#ylo.density                              = 1.225               
#ylo.temperature                          = 290.0               
#ylo.tke                                  = 0.0
#yhi.type                                 = pressure_outflow     

zlo.type                                 = wall_model
zhi.type                                 = slip_wall
zhi.temperature_type                     = fixed_gradient
zhi.temperature                          = 0.003

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
ICNS.source_terms                        = BoussinesqBuoyancy CoriolisForcing ABLForcing
##--------- Additions by calc_inflow_stats.py ---------#
#ABL.wall_shear_stress_type = "local"
#ABL.inflow_outflow_mode = true
#ABL.wf_velocity = XX XX
#ABL.wf_vmag = XX
#ABL.wf_theta = XX
#BodyForce.magnitude = XX XX XX
#BoussinesqBuoyancy.read_temperature_profile = true
#BoussinesqBuoyancy.tprofile_filename = avg_theta.dat
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
sampling.xy-domain.origin                = 0.0 0.0 86.5      
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
```
This file is almost identical to the spinup config file, except there are a few differences:
* I run this simulation with pre-defined regions where the mesh is refined. The refinement locations are specified in the `tagging` section, and I tell AMR-Wind to use `amr.max_level=2` levels of refinement. This means my finest grid cell is now 20 / 2 / 2 = 5 m wide.
* Because my finest mesh is now smaller, I also reduced my timestep by a factor of 4 to 0.125 seconds. 
* I am starting this simulation from the last timestep of the spinup simulation, using the `io.restart_file` line
* I am also now saving out boundary condition data, using the `ABL.bndry*` lines
* I plan to make some videos of hub-height planes and cross sections, so I have increased the frequency of `sampling.output_frequency`

Just like before, I create a projects directory `/projects/wakedynamics/orybchuk/amr-wind-tutorial/02_atmosphere/precursor` and a corresponding scratch directory `/scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/precursor`  . After linking the config file, I kick off a job with `sbatch runfile.sbatch`, and I confirm that this job successfully ran to completion by checking the log file.

Just like before, you can visualize the volume files using Paraview. Here, we're going to visualize the `sampling.xy-domain` and `sampling.xz-domain` data as well. That data is saved to `/scratch/orybchuk/wakedynamics/amr-wind-tutorial/02_atmosphere/precursor/post_processing/sampling14400.nc`. You can analyze data directly out of that file, but I like to reformat the data so that it is spatially sorted. To reformat data, run the code in `reformat_precursor_planes.ipynb`. To visualize the reformatted data, run `viz_precursor_planes.ipynb`.