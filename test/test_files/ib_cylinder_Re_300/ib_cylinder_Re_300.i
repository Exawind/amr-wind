#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   -10.0     # Max (simulated) time to evolve
time.max_step                =   20          # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
time.fixed_dt         =   -0.05        # Use this constant dt if > 0
time.cfl              =   0.45         # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval            =  10       # Steps between plot files
time.checkpoint_interval      =  -1       # Steps between checkpoint files

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
ConstValue.density.value = 1.0
ConstValue.velocity.value = 1.0 0.0 0.0

io.output_default_variables = 0
io.outputs = density p
io.derived_outputs = "components(velocity,0,1)" "components(gp,0,1)"

incflo.use_godunov = 1
incflo.godunov_type = "weno_z"
incflo.do_initial_proj = 1
incflo.initial_iterations = 3
transport.viscosity = 1.0e-3
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = Laminar

incflo.physics = FreeStream IB
IB.labels = IB1  
IB.IB1.type = Cylinder 
IB.IB1.center = 0.0 0.0 0.0
IB.IB1.radius = 0.05 
IB.IB1.height = 0.25

amr.n_cell     = 64 64 16   # Grid cells at coarsest AMRlevel
tagging.labels = sr                                                                                                                
tagging.sr.type = CartBoxRefinement                                                                                                                
tagging.sr.static_refinement_def = static_box.refine                                                                                             
amr.max_level = 2

geometry.prob_lo        =   -0.5 -0.5 -0.125
geometry.prob_hi        =    0.5  0.5  0.125  
geometry.is_periodic    =   0   0   1   # Periodicity x y z (0/1)

# Boundary conditions
xlo.type = "mass_inflow"
xlo.density = 1.0
xlo.velocity = 1.0 0.0 0.0
xhi.type = "pressure_outflow"
ylo.type =   "slip_wall"
yhi.type =   "slip_wall"

incflo.verbose          =   0          # incflo_level
nodal_proj.verbose = 0

nodal_proj.mg_rtol = 1.0e-10
nodal_proj.mg_atol = 1.0e-12
mac_proj.mg_rtol = 1.0e-10
mac_proj.mg_atol = 1.0e-12
