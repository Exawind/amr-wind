target_sources(${amr_wind_unit_test_exe_name} PRIVATE

  test_refinement.cpp
  test_plane_averaging.cpp
  test_field_plane_averaging.cpp
  test_second_moment.cpp
  test_sampling.cpp
  test_linear_interpolation.cpp
  test_free_surface.cpp
  test_wave_energy.cpp
  test_diagnostics.cpp
  )

if (AMR_WIND_ENABLE_NETCDF)
  target_sources(${amr_wind_unit_test_exe_name} PRIVATE
    test_ncutils.cpp
    )
endif()
