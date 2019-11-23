#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
stop_time               =   100          # Max (simulated) time to evolve
max_step                =   -1          # Max number of time steps
steady_state            =   0           # Steady-state solver? 

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
incflo.fixed_dt         =   -1          # Use this constant dt if > 0
incflo.cfl              =   0.49         # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
amr.plot_int            =   -1          # Steps between plot files
amr.plot_per            =   0.5         # Steps between plot files
amr.check_int           =   1000        # Steps between checkpoint files
amr.restart             =   ""          # Checkpoint to restart from 
amr.plt_p               =   1
amr.plt_tracer = 1

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity          =   0.  0.  -1.0  # Gravitational force (3D)
incflo.ro_0             =   1.          # Reference density 

incflo.fluid_model      =   "newtonian" # Fluid model (rheology)
incflo.mu               =   0.003       # Dynamic viscosity coefficient
incflo.mu_s             =   0.003
incflo.use_boussinesq   = 1
incflo.ntrac = 1
incflo.advect_tracer = 1

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              =   32 32  16  # Grid cells at coarsest AMRlevel
amr.max_level           =   0           # Max AMR level in hierarchy 

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   -1.0 -1.0 -0.5  # Lo corner coordinates
geometry.prob_hi        =   1.0 1.0 0.5 # Hi corner coordinates
geometry.is_periodic    =   1   1    0   # Periodicity x y z (0/1)

# Boundary conditions
zlo.type = "nsw"
zlo.velocity = 0. 0. 0.
zlo.tracer = 0.5

zhi.type = "nsw"
zhi.velocity = 0. 0. 0.
zhi.tracer = -0.5

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#           INITIAL CONDITIONS          #
#.......................................#
incflo.probtype         =   36
incflo.ic_u             =   0.0         #
incflo.ic_v             =   0.0         #
incflo.ic_w             =   0.01         #
incflo.ic_p             =   0.0         #

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          NUMERICAL PARAMETERS         #
#.......................................#
incflo.steady_state_tol =   1.e-4       # Tolerance for steady-state
amrex.fpe_trap_invalid  =   0           # Trap NaNs

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              VERBOSITY                #
#.......................................#
incflo.verbose          =   2           # incflo_level
mac.verbose             =   0           # MacProjector
