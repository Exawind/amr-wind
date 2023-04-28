target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      BCInterface.cpp
      FixedGradientBC.cpp
   )

add_subdirectory(wall_models)
