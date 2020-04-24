target_sources(${amr_wind_lib_name}
  PRIVATE

  ABL.cpp
  ABLFieldInit.cpp
  ABLForcing.cpp
  ABLWallFunction.cpp
  BoussinesqBuoyancy.cpp
  CoriolisForcing.cpp
  BoussinesqBubble.cpp
  BoussinesqBubbleFieldInit.cpp
  DensityBuoyancy.cpp
  RayleighTaylor.cpp
  RayleighTaylorFieldInit.cpp
  TaylorGreenVortex.cpp
  )

target_include_directories(
  ${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
