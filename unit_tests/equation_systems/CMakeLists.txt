target_sources(${amr_wind_unit_test_exe_name}
  PRIVATE

  test_pde.cpp
  )

target_include_directories(${amr_wind_unit_test_exe_name} PRIVATE
  ${CMAKE_SOURCE_DIR}/src/equation_systems)
