target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      boundary_conditions.cpp
      incflo_fillpatch.cpp
      incflo_set_density_bcs.cpp
      incflo_set_tracer_bcs.cpp
      incflo_set_velocity_bcs.cpp
      #Headers
      bc_mod_F.H
      #Fortran
      bc_mod.f90
      fill_bc0.f90
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
