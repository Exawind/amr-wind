target_sources(${amr_wind_unit_test_exe_name} PRIVATE
  # testing utilities
  abl_test_utils.cpp

  # test cases
  test_abl_init.cpp
  test_abl_src.cpp
  )

if (AMR_WIND_ENABLE_NETCDF)
  target_sources(${amr_wind_unit_test_exe_name} PRIVATE
    test_abl_init_ncf.cpp
    )
endif()

add_subdirectory(actuator)
