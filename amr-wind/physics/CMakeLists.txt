target_sources(${amr_wind_lib_name}
  PRIVATE

  BoussinesqBubble.cpp
  BoussinesqBubbleFieldInit.cpp
  RayleighTaylor.cpp
  RayleighTaylorFieldInit.cpp
  TaylorGreenVortex.cpp
  FreeStream.cpp
  ConvectingTaylorVortex.cpp
  )
