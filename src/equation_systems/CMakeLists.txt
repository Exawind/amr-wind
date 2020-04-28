target_sources(${amr_wind_lib_name} PRIVATE
  PDEBase.cpp
  DiffusionOps.cpp
  )

add_subdirectory(icns)
add_subdirectory(temperature)
add_subdirectory(density)
add_subdirectory(tke)
add_subdirectory(sdr)

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
