target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      diagnostics.cpp
      incflo_build_info.cpp
      io.cpp
      PlaneAveraging.cpp
      bc_ops.cpp
   )

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(tagging)
