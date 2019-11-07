#Aggregate AMR-Wind source files
function(get_amr_wind_sources amr_wind_exe_name)
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/advance.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_compute_dt.cpp
     ${AMR_WIND_SOURCE_DIR}/main.cpp
     ${AMR_WIND_SOURCE_DIR}/param_mod.f90
  )
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/boundary_conditions")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/bc_mod.f90
     ${AMR_WIND_SOURCE_DIR}/fill_bc0.f90
     ${AMR_WIND_SOURCE_DIR}/boundary_conditions.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_fillpatch.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_set_density_bcs.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_set_tracer_bcs.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_set_velocity_bcs.cpp
  )
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/convection")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/incflo_compute_convective_term.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_compute_fluxes.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_MAC_projection.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_predict_vels_on_faces.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_set_mac_velocity_bcs.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_slopes.cpp
  )
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/derive")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/derive_mod.f90
     ${AMR_WIND_SOURCE_DIR}/derive.cpp
  )
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/diffusion")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/DiffusionOp.cpp
     ${AMR_WIND_SOURCE_DIR}/set_scal_diff_bcs.f90
     ${AMR_WIND_SOURCE_DIR}/set_vel_diff_bcs.f90
  )
  if(AMR_WIND_ENABLE_EB)
     set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/embedded_boundaries")
     add_sources(GlobalSourceList
        ${AMR_WIND_SOURCE_DIR}/get_eb_walls.f90
        ${AMR_WIND_SOURCE_DIR}/embedded_boundaries.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_annulus.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_box.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_cylinder.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_regular.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_sphere.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_spherecube.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_tuscan.cpp
        ${AMR_WIND_SOURCE_DIR}/eb_twocylinders.cpp
        ${AMR_WIND_SOURCE_DIR}/get_walls.cpp
        ${AMR_WIND_SOURCE_DIR}/writeEBsurface.cpp
     )
  endif()
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/projection")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/incflo_apply_nodal_projection.cpp
  )
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/rheology")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/rheology_mod.f90
     ${AMR_WIND_SOURCE_DIR}/rheology.cpp
  )
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/setup")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/incflo_to_fortran.F90
     ${AMR_WIND_SOURCE_DIR}/init_fluid.f90 
     ${AMR_WIND_SOURCE_DIR}/set_bc_type.f90 
     ${AMR_WIND_SOURCE_DIR}/set_delp_dir.f90
     ${AMR_WIND_SOURCE_DIR}/set_p0.f90 
     ${AMR_WIND_SOURCE_DIR}/set_ppe_bcs.f90 
     ${AMR_WIND_SOURCE_DIR}/incflo_arrays.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_setup_solvers.cpp
     ${AMR_WIND_SOURCE_DIR}/init.cpp
  )
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/utilities")
  add_sources(GlobalSourceList
     ${AMR_WIND_SOURCE_DIR}/constant_mod.f90
     ${AMR_WIND_SOURCE_DIR}/diagnostics.cpp
     ${AMR_WIND_SOURCE_DIR}/incflo_build_info.cpp
     ${AMR_WIND_SOURCE_DIR}/io.cpp
  )
  #Add generated source files
  add_sources(GlobalSourceList
     #${CMAKE_BINARY_DIR}/generated_files/${amr_wind_exe_name}_generated_files/extern.f90
     ${CMAKE_BINARY_DIR}/generated_files/${amr_wind_exe_name}_generated_files/AMReX_buildInfo.cpp
  )
endfunction(get_amr_wind_sources amr_wind_exe_name)
