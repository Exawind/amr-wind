target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      rheology.cpp
      #Headers
      rheology_F.H
      #Fortran
      rheology_mod.f90
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
