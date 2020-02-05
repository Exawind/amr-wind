target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      boundary_conditions.cpp
      incflo_fillpatch.cpp
      incflo_fillphysbc.cpp
      incflo_set_velocity_bcs.cpp
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
