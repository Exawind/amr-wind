target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      incflo_arrays.cpp
      init.cpp
      set_background_pressure.cpp      
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
