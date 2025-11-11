target_sources(${amr_wind_unit_test_exe_name}
  PRIVATE
  utest_main.cpp
  test_config.cpp
  )

add_subdirectory(aw_test_utils)
add_subdirectory(core)
add_subdirectory(wind_energy)
add_subdirectory(offshore_wind)
add_subdirectory(utilities)
add_subdirectory(equation_systems)
add_subdirectory(turbulence)
add_subdirectory(fvm)
add_subdirectory(multiphase)
add_subdirectory(ocean_waves)
add_subdirectory(projection)
add_subdirectory(boundary_conditions)
add_subdirectory(eos_models)

if(AMR_WIND_ENABLE_MASA)
  add_subdirectory(mms)
endif()

target_compile_options(${amr_wind_unit_test_exe_name} PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:${AMR_WIND_CXX_FLAGS}>)
target_include_directories(${amr_wind_unit_test_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#target_link_libraries(${utest_exe_name} PRIVATE GTest::GTest GTest::Main)
target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE gtest)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE
  ${PROJECT_SOURCE_DIR}/submods/googletest/googletest/include)

target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE ${amr_wind_lib_name} AMReX-Hydro::amrex_hydro_api)
if (AMR_WIND_ENABLE_W2A)
  target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE Waves2AMR::Waves2AMR)
endif()

install(TARGETS ${amr_wind_unit_test_exe_name})
