target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      embedded_boundaries.cpp
      eb_cyl_tuscan.cpp
      eb_annulus.cpp
      eb_box.cpp
      eb_cylinder.cpp
      eb_regular.cpp
      eb_sphere.cpp
      eb_spherecube.cpp
      eb_tuscan.cpp
      eb_twocylinders.cpp
      writeEBsurface.cpp
      #Headers
      eb_if.H
      #Fortran
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
