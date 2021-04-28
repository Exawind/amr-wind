time.stop_time               =   6     # Max (simulated) time to evolve
time.max_step                =   20          # Max number of time steps

time.fixed_dt         =   0.001        # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor

time.plot_interval            =  10       # Steps between plot files
time.checkpoint_interval      =  -100       # Steps between checkpoint files

io.output_default_variables = 0
io.outputs = density vof
io.derived_outputs = "components(velocity,0,2)" "components(gp,0,2)"

incflo.use_godunov = 1
incflo.godunov_type = "weno"
transport.model = TwoPhaseTransport
transport.viscosity_fluid1=1.e-6
transport.viscosity_fluid2=1.48e-5

transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Laminar 

incflo.physics = MultiPhase DamBreak 
MultiPhase.density_fluid1=1000.
MultiPhase.density_fluid2=1.
ICNS.source_terms = GravityForcing 

amr.n_cell              = 64 16 64    # Grid cells at coarsest AMRlevel
amr.max_level           = 0           # Max AMR level in hierarchy 
amr.blocking_factor     = 8

geometry.prob_lo        =   0   0.   0.  # Lo corner coordinates
geometry.prob_hi        =   0.5   0.125  0.5  # Hi corner coordinates
geometry.is_periodic    =   0   1   0   # Periodicity x y z (0/1)

xlo.type =   "slip_wall"
xhi.type =   "slip_wall"
zlo.type =   "slip_wall"
zhi.type =   "slip_wall"

incflo.verbose=0
