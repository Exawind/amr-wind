target_sources(${amr_wind_lib_name}
  PRIVATE
  MultiPhase.cpp
  VortexPatch.cpp
  ZalesakDisk.cpp
  DamBreak.cpp
  SloshingTank.cpp
  RainDrop.cpp
  )

add_subdirectory(wave_basin)

