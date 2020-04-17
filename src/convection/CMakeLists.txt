target_sources(${amr_wind_lib_name}
   PRIVATE
      #C++
      incflo_MAC_bcs.cpp
      incflo_godunov_ppm.cpp
      incflo_MAC_projection.cpp	
      incflo_godunov_plm.cpp	
      incflo_godunov_predict.cpp
      incflo_mol_predict.cpp
      incflo_godunov_advection.cpp
      incflo_mol_fluxes.cpp
   )

target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
