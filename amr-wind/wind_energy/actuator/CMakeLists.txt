target_sources(${amr_wind_lib_name}
  PRIVATE

  actuator_utils.cpp
  Actuator.cpp
  ActuatorContainer.cpp
  )

add_subdirectory(wing)
