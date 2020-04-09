target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      init.cpp
      set_background_pressure.cpp      
   )

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
