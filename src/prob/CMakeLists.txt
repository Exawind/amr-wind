target_sources(${amr_wind_lib_name}
  PRIVATE
  prob_bc.cpp
  prob_init_fluid.cpp
  )

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
