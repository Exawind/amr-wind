target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      diagnostics.cpp
      incflo_build_info.cpp
      incflo_steady_state.cpp
      io.cpp
      PlaneAveraging.cpp
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
