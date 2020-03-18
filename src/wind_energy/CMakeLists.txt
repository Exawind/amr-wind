target_sources(${amr_wind_exe_name}
  PRIVATE

  ABLFieldInit.cpp
  ABL.cpp
  )

target_include_directories(
  ${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
