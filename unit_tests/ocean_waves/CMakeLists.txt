target_sources(${amr_wind_unit_test_exe_name} PRIVATE
  # test cases
  test_wave_utils.cpp
  test_relaxation_zones.cpp
  test_wave_theories.cpp
  )

if (AMR_WIND_ENABLE_W2A)
  target_sources(${amr_wind_unit_test_exe_name} PRIVATE
    test_waves_2_amr.cpp
    )
endif()