time.stop_time               = 0.5
time.fixed_dt               = .0017578125

incflo.use_godunov                     = true
incflo.godunov_type                     = ppm

time.plot_interval            = 9
time.checkpoint_interval           =  -100         # Steps between checkpoint files

transport.viscosity = 0.0
turbulence.model = Laminar

tagging.labels = static
tagging.static.type = CartBoxRefinement
tagging.static.static_refinement_def = static_box.txt

amr.max_level           =   1
amr.n_cell              = 128 4 4
amr.max_grid_size       = 4
amr.blocking_factor     = 4
amr.n_error_buf     =   0            # Blocking factor for grids

geometry.prob_lo        = 0.0 -0.0156250000 -0.0156250000
geometry.prob_hi        = 1.0 0.0156250000 0.0156250000

geometry.is_periodic    =   1   1   1   # Periodicity x y z (0/1)

incflo.physics = ScalarAdvection
scalaradvection.u = 1
scalaradvection.v = 0
scalaradvection.x0 = 0.25
scalaradvection.y0 = 0.5
scalaradvection.amplitude = 1.0
scalaradvection.x_width = 0.0251646060
scalaradvection.y_width = 0
scalaradvection.x_wavenumber = 226.1946710304
scalaradvection.y_wavenumber = 0
scalaradvection.shape = gaussianwavepacket
scalaradvection.output_fname = error.log

incflo.diffusion_type   = 2             # 0 = Explicit, 1 = Crank-Nicolson, 2 = Implicit

incflo.verbose          =   1           # incflo_level
mac_proj.verbose        =   0           # MAC Projector
nodal_proj.verbose      =   0           # Nodal Projector
time.use_force_cfl = false
