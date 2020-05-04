target_sources(${amr_wind_lib_name}
  PRIVATE

  BoussinesqBubble.cpp
  BoussinesqBubbleFieldInit.cpp
  RayleighTaylor.cpp
  RayleighTaylorFieldInit.cpp
  TaylorGreenVortex.cpp
  FreeStream.cpp
  )

target_include_directories(
  ${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
