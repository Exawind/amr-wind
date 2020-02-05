target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      incflo_compute_advection_term.cpp
      incflo_compute_fluxes.cpp
      incflo_MAC_projection.cpp
      incflo_predict_vels_on_faces.cpp
      incflo_convection_K.cpp
      incflo_MAC_bcs.cpp
      incflo_godunov_predict.cpp
      incflo_godunov_advection.cpp
      #Headers
      incflo_godunov_ppm.H
      incflo_MAC_bcs.H
      incflo_convection_K.H
   )

target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
