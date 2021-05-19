target_sources(${amr_wind_lib_name}
  PRIVATE

  ABL.cpp
  ABLStats.cpp
  ABLFieldInit.cpp
  ABLWallFunction.cpp
  ABLFillInflow.cpp
  ABLBoundaryPlane.cpp
  MOData.cpp
  ABLWrf.cpp
  ABLStatsFileRead.cpp
  )

add_subdirectory(actuator)
