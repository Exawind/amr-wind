target_sources(${amr_wind_lib_name} PRIVATE
  ABLForcing.cpp
  BoussinesqBuoyancy.cpp
  DensityBuoyancy.cpp
  CoriolisForcing.cpp
  )

target_include_directories(
${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
