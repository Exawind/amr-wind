target_sources(${amr_wind_lib_name}
  PRIVATE

  ABL.cpp
  ABLStats.cpp
  ABLFieldInit.cpp
  ABLFieldInitFile.cpp
  ABLWallFunction.cpp
  ABLFillInflow.cpp
  ABLBoundaryPlane.cpp
  MOData.cpp
  ABLFillMPL.cpp
  ABLModulatedPowerLaw.cpp
  )

add_subdirectory(actuator)
