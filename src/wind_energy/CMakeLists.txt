target_sources(${amr_wind_exe_name}
  PRIVATE

  ABL.cpp
  ABLFieldInit.cpp
  ABLForcing.cpp
  ABLWallFunction.cpp
  BoussinesqBuoyancy.cpp
  CoriolisForcing.cpp
  )

target_include_directories(
  ${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
