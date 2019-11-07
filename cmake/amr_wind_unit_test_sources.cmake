#Aggregate AMR-Wind unit test source files
function(get_amr_wind_unit_test_sources)
  set(AMR_WIND_SOURCE_DIR "${CMAKE_SOURCE_DIR}/test/test_files/unit-tests")
  add_sources(GlobalUnitSourceList
     ${AMR_WIND_SOURCE_DIR}/unit-tests-${AMR_WIND_DIM}d.C
     ${AMR_WIND_SOURCE_DIR}/unit-test-${AMR_WIND_DIM}d-1.C
  )
endfunction(get_amr_wind_unit_test_sources)
