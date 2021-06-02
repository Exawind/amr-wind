target_sources(${amr_wind_lib_name} PRIVATE
  FastIface.cpp
  )

if (AMR_WIND_ENABLE_OPENFAST)
  target_sources(${amr_wind_lib_name} PRIVATE
    TurbineFast.cpp
    )
endif()
