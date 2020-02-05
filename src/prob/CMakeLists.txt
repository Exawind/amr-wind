target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      prob_bc.cpp 
      prob_init_fluid.cpp
      #Headers
      prob_bc.H
      incflo_prob_I.H      
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
