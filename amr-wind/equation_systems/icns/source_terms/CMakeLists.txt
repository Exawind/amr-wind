target_sources(${amr_wind_lib_name} PRIVATE
  ABLForcing.cpp
  GeostrophicForcing.cpp
  IdealizedHurricaneForcing.cpp
  BoussinesqBuoyancy.cpp
  DensityBuoyancy.cpp
  GravityForcing.cpp
  CoriolisForcing.cpp
  BodyForce.cpp
  ABLMeanBoussinesq.cpp
  ActuatorForcing.cpp
  SynthTurbForcing.cpp
  )
