target_sources(${amr_wind_lib_name}
  PRIVATE
  incflo_advance.cpp
  incflo.cpp
  incflo_compute_dt.cpp
  incflo_tagging.cpp
  incflo_regrid.cpp
  incflo_compute_forces.cpp
  incflo_field_repo.cpp
  )

target_sources(${amr_wind_exe_name} PRIVATE main.cpp)

#AMR-Wind include directories
target_include_directories(${amr_wind_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

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
add_subdirectory(wind_energy)
add_subdirectory(equation_systems)
add_subdirectory(transport_models)
add_subdirectory(turbulence)

include(AMReXBuildInfo)
generate_buildinfo(${amr_wind_lib_name} ${AMREX_SUBMOD_LOCATION})
target_include_directories(${amr_wind_lib_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

target_link_libraries_system(${amr_wind_lib_name} PUBLIC amrex)
target_link_libraries(${amr_wind_exe_name} PRIVATE ${amr_wind_lib_name})
