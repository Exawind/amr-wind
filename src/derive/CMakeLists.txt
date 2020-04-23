target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      derive.cpp
      compute_strainrate.cpp
   )

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
