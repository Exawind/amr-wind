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
  ABLMesoscaleForcing.cpp
  ABLMesoscaleInput.cpp
  ABLFillMPL.cpp
  ABLModulatedPowerLaw.cpp
  ABLAnelastic.cpp
  )

add_subdirectory(actuator)
