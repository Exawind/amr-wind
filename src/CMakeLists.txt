target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      main.cpp
      incflo_advance.cpp
      incflo.cpp
      incflo_compute_dt.cpp
      incflo_tagging.cpp
      incflo_regrid.cpp
      incflo_compute_forces.cpp
      #Headers
   )

#Add generated source files
set_property(SOURCE ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp PROPERTY GENERATED 1)
target_sources(${amr_wind_exe_name}
   PRIVATE
      ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
)

#AMR-Wind include directories
target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#Needed for AMReX_buildInfo.H
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Tools/C_scripts)

#Gather all other source files
add_subdirectory(core)
add_subdirectory(boundary_conditions)
add_subdirectory(convection)
add_subdirectory(derive)
add_subdirectory(diffusion)
add_subdirectory(projection)
add_subdirectory(rheology)
add_subdirectory(setup)
add_subdirectory(utilities)
add_subdirectory(prob)
if(AMR_WIND_ENABLE_EB)
   add_subdirectory(embedded_boundaries)
endif()

#Link to amrex library
target_link_libraries_system(${amr_wind_exe_name} amrex)
