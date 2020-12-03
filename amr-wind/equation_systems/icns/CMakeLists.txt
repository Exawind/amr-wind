target_sources(${amr_wind_lib_name} PRIVATE
  icns_ops.cpp
  icns_advection.cpp
  icns.cpp)

add_subdirectory(source_terms)
