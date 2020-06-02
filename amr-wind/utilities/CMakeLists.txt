target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      diagnostics.cpp
      io.cpp
      PlaneAveraging.cpp
      bc_ops.cpp
      console_io.cpp
      IOManager.cpp
      FieldPlaneAveraging.cpp
      SecondMomentAveraging.cpp
      ThirdMomentAveraging.cpp
   )

add_subdirectory(tagging)
