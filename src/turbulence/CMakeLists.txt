target_sources(${amr_wind_lib_name} PRIVATE
  LaminarModel.cpp
  )

add_subdirectory(LES)

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
