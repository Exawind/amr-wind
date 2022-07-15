target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      diagnostics.cpp
      io.cpp
      bc_ops.cpp
      console_io.cpp
      IOManager.cpp
      FieldPlaneAveraging.cpp
      SecondMomentAveraging.cpp
      ThirdMomentAveraging.cpp

      PostProcessing.cpp
      DerivedQuantity.cpp
      DerivedQtyDefs.cpp
   )

add_subdirectory(tagging)
add_subdirectory(sampling)
add_subdirectory(averaging)

if (AMR_WIND_ENABLE_NETCDF)
  add_subdirectory(ncutils)
endif()

if (AMR_WIND_ENABLE_ASCENT)
   add_subdirectory(ascent)
endif()
