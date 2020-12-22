target_sources(${amr_wind_unit_test_exe_name} PRIVATE
  # testing utilities
  abl_test_utils.cpp

  # test cases
  test_abl_init.cpp
  test_abl_src.cpp
  )

add_subdirectory(actuator)
