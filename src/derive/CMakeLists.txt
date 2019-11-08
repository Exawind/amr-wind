target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      derive.cpp
      #Headers
      derive_F.H
      #Fortran
      derive_mod.f90
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
