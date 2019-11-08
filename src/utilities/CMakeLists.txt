target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      diagnostics.cpp
      incflo_build_info.cpp
      io.cpp
      #Fortran
      constant_mod.f90
   )

#target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
