time.stop_time               =   20     # Max (simulated) time to evolve
time.max_step                =   10          # Max number of time steps

time.fixed_dt         =   -0.005        # Use this constant dt if > 0
time.cfl              =   0.5         # CFL factor
time.initial_dt =0.001
time.plot_interval            =  10       # Steps between plot files
time.checkpoint_interval      =  10       # Steps between checkpoint files

io.output_default_variables = 0
io.outputs = density vof
io.derived_outputs = "components(velocity,0,2)" "components(gp,0,2)"

incflo.use_godunov = 1
incflo.godunov_type="weno"
transport.model = TwoPhaseTransport
transport.viscosity_fluid1=0.03132
transport.viscosity_fluid2=0.000018
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Laminar 

incflo.physics = MultiPhase SloshingTank 
SloshingTank.amplitude=0.05
MultiPhase.density_fluid1=1000.
MultiPhase.density_fluid2=1.
ICNS.source_terms = GravityForcing 

amr.n_cell              = 128 16 96    # Grid cells at coarsest AMRlevel
tagging.labels = sr
tagging.sr.type = CartBoxRefinement
tagging.sr.static_refinement_def = static_box.refine
amr.max_level = 1 

geometry.prob_lo        =   -1   0.   -1.   # Lo corner coordinates
geometry.prob_hi        =    1   0.25  0.5  # Hi corner coordinates
geometry.is_periodic    =    0   1    0   # Periodicity x y z (0/1)

xlo.type =   "slip_wall"
xhi.type =   "slip_wall"
zlo.type =   "slip_wall"
zhi.type =   "slip_wall"

incflo.verbose          =   0          # incflo_level

