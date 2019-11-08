target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      incflo_arrays.cpp
      incflo_setup_solvers.cpp
      init.cpp
      #Headers
      setup_F.H
      #Fortran
      incflo_to_fortran.F90
      init_fluid.f90 
      set_bc_type.f90 
      set_delp_dir.f90
      set_p0.f90 
      set_ppe_bcs.f90 
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
