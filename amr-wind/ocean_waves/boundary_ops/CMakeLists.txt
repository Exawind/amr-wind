target_sources(${amr_wind_lib_name}
  PRIVATE
  OceanWaves.cpp
  OceanWavesFillInflow.cpp
  )
add_subdirectory(relaxation_zones)
