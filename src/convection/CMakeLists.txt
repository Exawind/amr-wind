target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      incflo_compute_convective_term.cpp
      incflo_compute_fluxes.cpp
      incflo_MAC_projection.cpp
      incflo_predict_vels_on_faces.cpp
      incflo_set_mac_velocity_bcs.cpp
      incflo_slopes.cpp
   )

#target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
