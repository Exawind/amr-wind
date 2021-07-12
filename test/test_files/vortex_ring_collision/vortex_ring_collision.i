#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   22.0   # Max (simulated) time to evolve
time.max_step                =   -20   # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
#time.fixed_dt         =   0.006835938        # Use this constant dt if > 0
time.cfl              =   0.5        # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval  =  100   # Steps between plot files
time.checkpoint_interval =   -1  # Steps between checkpoint files
io.output_default_variables = 0
io.outputs = density p
io.derived_outputs = q_criterion q_criterion_nondim mag_vorticity
incflo.post_processing = ke enst
ke.type = KineticEnergy
enst.type = Enstrophy

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.use_godunov = 1
incflo.godunov_type = "plm"
transport.viscosity = 0.001
transport.laminar_prandtl = 1.0
turbulence.model = Laminar


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              =  400 400 200   # Grid cells at coarsest AMRlevel
#amr.max_grid_size       =   64 64 64   
amr.max_level           =   4           # Max AMR level in hierarchy 
#amr.blocking_factor     =   50 50 50

tagging.labels = qc
tagging.qc.type = VorticityMagRefinement
tagging.qc.nondim = true
tagging.qc.values = 0.1 0.1 0.1 0.1
time.regrid_interval = 10

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   -20. -20. -10.  # Lo corner coordinates
geometry.prob_hi        =   20.  20.  10.  # Hi corner coordinates
geometry.is_periodic    =   0   0   1   # Periodicity x y z (0/1)
xlo.type="pressure_outflow"
ylo.type="pressure_outflow"
xhi.type="pressure_outflow"
yhi.type="pressure_outflow"

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          INITIAL CONDITIONS           #
#.......................................#
incflo.physics = VortexRingCollision
vortexringcollision.R = 1.0
vortexringcollision.Gamma = 1.0
vortexringcollision.delta = 0.1
vortexringcollision.dz = 2.0
