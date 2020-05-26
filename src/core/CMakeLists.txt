target_sources(${amr_wind_lib_name}
  PRIVATE
  SimTime.cpp

  Field.cpp
  IntField.cpp
  FieldRepo.cpp
  )

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
