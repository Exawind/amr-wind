#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
stop_time               =   10          # Max (simulated) time to evolve
max_step                =   -1          # Max number of time steps
steady_state            =   0           # Steady-state solver? 

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
incflo.fixed_dt         =   -1          # Use this constant dt if > 0
incflo.cfl              =   0.9         # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
amr.plot_int            =   -1          # Steps between plot files
amr.plot_per            =   0.5         # Steps between plot files
amr.check_int           =   1000        # Steps between checkpoint files
amr.restart             =   ""          # Checkpoint to restart from 
amr.plt_p               =   1

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.gravity          =   0.  0.  0.  # Gravitational force (3D)
incflo.ro_0             =   1.          # Reference density 

incflo.fluid_model      =   "newtonian" # Fluid model (rheology)
incflo.mu               =   0.001       # Dynamic viscosity coefficient

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              =   128 32  16  # Grid cells at coarsest AMRlevel
amr.max_level           =   0           # Max AMR level in hierarchy 

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   0.  0.   0.  # Lo corner coordinates
geometry.prob_hi        =   1.64 .41 .205 # Hi corner coordinates
geometry.is_periodic    =   0   0    1   # Periodicity x y z (0/1)

# Add cylinder 
incflo.geometry         =   "cylinder"
cylinder.internal_flow  =   false
cylinder.radius         =   0.05
cylinder.direction      =   2
cylinder.center         =   .2  .2  0. 

# Boundary conditions
xlo.type                =   "mi"
xlo.velocity            =   0.2 0.  0. 
xhi.type                =   "po"
xhi.pressure            =   0.0
ylo.type                =   "nsw"
ylo.velocity            =   0.  0.  0.
yhi.type                =   "nsw"
yhi.velocity            =   0.  0.  0.

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#           INITIAL CONDITIONS          #
#.......................................#
incflo.probtype         =   3
incflo.ic_u             =   0.2         #
incflo.ic_v             =   0.0         #
incflo.ic_w             =   0.0         #
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
