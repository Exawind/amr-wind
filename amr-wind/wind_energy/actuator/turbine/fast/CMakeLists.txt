if (AMR_WIND_ENABLE_OPENFAST)
  find_package(OpenFAST REQUIRED)

  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_OPENFAST)
  target_include_directories(${amr_wind_lib_name} PUBLIC ${OpenFAST_INCLUDE_DIRS})
  target_link_libraries(${amr_wind_lib_name} PUBLIC ${OpenFAST_LIBRARIES})
endif()

target_sources(${amr_wind_lib_name} PRIVATE
  FastIface.cpp
  )

if (AMR_WIND_ENABLE_OPENFAST)
  target_sources(${amr_wind_lib_name} PRIVATE
    TurbineFast.cpp
    )
endif()
