target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      incflo_diffusion.cpp
      DiffusionTensorOp.cpp
      DiffusionScalarOp.cpp
      #Headers
      DiffusionTensorOp.H 
      DiffusionScalarOp.H
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
