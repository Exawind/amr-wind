target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      diagnostics.cpp
      incflo_build_info.cpp
      io.cpp
      PlaneAveraging.cpp
      bc_ops.cpp
      console_io.cpp
      IOManager.cpp
   )

add_subdirectory(tagging)
