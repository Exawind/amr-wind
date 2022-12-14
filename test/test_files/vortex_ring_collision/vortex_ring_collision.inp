#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
#time.stop_time               =   22.0   # Max (simulated) time to evolve
time.max_step                =   20   # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
#time.fixed_dt         =   0.006835938        # Use this constant dt if > 0
time.cfl              =   0.5        # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval  =  20   # Steps between plot files
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
amr.n_cell              =  32 32 16   # Grid cells at coarsest AMRlevel
#amr.max_grid_size       =   64 64 64   
amr.max_level           =   1           # Max AMR level in hierarchy 
amr.blocking_factor     =   4 4 4

tagging.labels = vm geom

tagging.vm.type = VorticityMagRefinement
tagging.vm.nondim = true
tagging.vm.values = 0.1 0.1 0.1 0.1
time.regrid_interval = 10

tagging.geom.type = GeometryRefinement
tagging.geom.shapes = c1 c2
tagging.geom.min_level = 0
tagging.geom.max_level = 10
tagging.geom.c1.type = cylinder
tagging.geom.c1.start = 0.0 0.0 -1.5
tagging.geom.c1.end = 0.0 0.0 -0.5
tagging.geom.c1.inner_radius = 0.5
tagging.geom.c1.outer_radius = 1.5

tagging.geom.c2.type = cylinder
tagging.geom.c2.start = 0.0 0.0 0.5
tagging.geom.c2.end = 0.0 0.0 1.5
tagging.geom.c2.inner_radius = 0.5
tagging.geom.c2.outer_radius = 1.5

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
incflo.physics = VortexRing
vortexring.type = collidingrings
vortexring.R = 1.0
vortexring.Gamma = 1.0
vortexring.delta = 0.1
vortexring.dz = 2.0
