#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
#time.stop_time               =   8.00   # Max (simulated) time to evolve
time.max_step                =   20   # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   0.04        # Use this constant dt if > 0
#time.cfl              =   0.5        # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval  =  32   # Steps between plot files
time.checkpoint_interval =   32  # Steps between checkpoint files
io.outputs = vorticity
io.derived_outputs = q_criterion q_criterion_nondim mag_vorticity
incflo.post_processing = ke enst
ke.type = KineticEnergy
ke.output_frequency = 1
enst.type = Enstrophy
enst.output_frequency = 1

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
amr.n_cell              =   32 32 32    # Grid cells at coarsest AMRlevel
#amr.max_grid_size       =   64 64 64   
amr.max_level           =   1           # Max AMR level in hierarchy 
amr.blocking_factor     =   16 16 16

tagging.labels = vm geom
tagging.vm.type = VorticityMagRefinement
tagging.vm.nondim = true
tagging.vm.values = 0.01 0.01 0.01 0.01
time.regrid_interval = 1

tagging.geom.type = GeometryRefinement
tagging.geom.shapes = c1
tagging.geom.min_level = 0
tagging.geom.max_level = 10
tagging.geom.c1.type = cylinder
tagging.geom.c1.start = 0.0 0.0 -0.5
tagging.geom.c1.end = 0.0 0.0 0.5
tagging.geom.c1.inner_radius = 0.5
tagging.geom.c1.outer_radius = 1.5

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   -6. -6. -6.  # Lo corner coordinates
geometry.prob_hi        =   6.  6.  6.  # Hi corner coordinates
geometry.is_periodic    =   1   1   1   # Periodicity x y z (0/1)

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          INITIAL CONDITIONS           #
#.......................................#
incflo.physics = VortexRing
vortexring.type = fatcore
vortexring.R = 1.0
vortexring.Gamma = 1.0
