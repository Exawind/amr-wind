target_sources(${amr_wind_lib_name} PRIVATE
  FastIface.cpp
  )

if (AMR_WIND_ENABLE_OPENFAST)
  add_subdirectory(fast_inp)
  target_sources(${amr_wind_lib_name} PRIVATE
    TurbineFast.cpp
    )
endif()
