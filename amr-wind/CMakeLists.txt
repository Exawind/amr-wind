target_sources(${amr_wind_lib_name}
  PRIVATE
  incflo_advance.cpp
  incflo.cpp
  incflo_compute_dt.cpp
  incflo_regrid.cpp
  helics.cpp
  CFDSim.cpp
  )

target_sources(${amr_wind_exe_name} PRIVATE main.cpp)

#AMR-Wind include directories
target_include_directories(${amr_wind_lib_name} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

#Gather all other source files
add_subdirectory(core)
add_subdirectory(boundary_conditions)
add_subdirectory(convection)
add_subdirectory(derive)
add_subdirectory(diffusion)
add_subdirectory(projection)
add_subdirectory(setup)
add_subdirectory(utilities)
add_subdirectory(wind_energy)
add_subdirectory(equation_systems)
add_subdirectory(transport_models)
add_subdirectory(turbulence)
add_subdirectory(physics)
add_subdirectory(overset)
add_subdirectory(immersed_boundary)
add_subdirectory(mesh_mapping_models)

include(AMReXBuildInfo)
generate_buildinfo(${amr_wind_lib_name} ${CMAKE_SOURCE_DIR})

# Generate AMR-Wind version header
configure_file("${CMAKE_SOURCE_DIR}/cmake/AMRWindVersion.H.in"
  "${CMAKE_CURRENT_BINARY_DIR}/AMRWindVersion.H" @ONLY)

target_link_libraries_system(${amr_wind_lib_name} PUBLIC AMReX::amrex AMReX-Hydro::amrex_hydro_api)
target_link_libraries(${amr_wind_exe_name} PRIVATE ${amr_wind_lib_name} AMReX-Hydro::amrex_hydro_api)
target_link_libraries(${aw_api_lib} PUBLIC ${amr_wind_lib_name} AMReX-Hydro::amrex_hydro_api)

# Set -fpic options
set_target_properties(${amr_wind_lib_name} buildInfo${amr_wind_lib_name}
  PROPERTIES POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE amr_wind_header_list ${CMAKE_CURRENT_SOURCE_DIR} "*.H")
add_custom_target(amr_wind_headers SOURCES ${amr_wind_header_list})

# Add CMake-generated headers to installation
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/AMRWindVersion.H
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
