#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   3     # Max (simulated) time to evolve
time.max_step                =   100          # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   -0.005        # Use this constant dt if > 0
time.cfl              =   0.45         # CFL factor
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval            =  50       # Steps between plot files
time.checkpoint_interval      =  100       # Steps between checkpoint files

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.use_godunov = 1
incflo.godunov_type="weno_z"
transport.model = TwoPhaseTransport
transport.viscosity_fluid1=0.0
transport.viscosity_fluid2=0.0
turbulence.model = KOmegaSST
TKE.source_terms = KwSSTSrc
SDR.source_terms = SDRSrc

incflo.physics = MultiPhase BreakingWaves 
BreakingWaves.amplitude=0.112
MultiPhase.density_fluid1=998.
MultiPhase.density_fluid2=1.2
ICNS.source_terms = GravityForcing 
MultiPhase.verbose=1

io.outputs = density velocity_mueff sdr tke mu_turb 
io.derived_outputs = "components(velocity,0)"
#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              = 64 16 48    # Grid cells at coarsest AMRlevel
amr.max_level = 0 

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =     0.0   0.0  -0.5   # Lo corner coordinates
geometry.prob_hi        =     2.0   0.5   1.0  # Hi corner coordinates
geometry.is_periodic    =     1     1     0   # Periodicity x y z (0/1)

zlo.type =   "slip_wall"
zhi.type =   "slip_wall"

incflo.verbose=1

