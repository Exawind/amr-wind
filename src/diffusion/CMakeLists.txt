target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      DiffusionOp.cpp
      #Headers
      DiffusionOp.H
      diffusion_F.H
      #Fortran
      set_scal_diff_bcs.f90
      set_vel_diff_bcs.f90
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
