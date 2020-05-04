target_sources(${amr_wind_lib_name}
  PRIVATE

  ABL.cpp
  ABLFieldInit.cpp
  ABLWallFunction.cpp
  )

target_include_directories(
  ${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
