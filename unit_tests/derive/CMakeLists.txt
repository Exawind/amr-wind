target_sources(
  ${amr_wind_unit_test_exe_name} PRIVATE

  test_strainrate.cpp
  AnalyticalFunctions.cpp
  test_scalar_gradient.cpp
  test_vector_gradient.cpp

  )
