target_sources(${amr_wind_lib_name}
  PRIVATE
  incflo_advance.cpp
  incflo.cpp
  incflo_compute_dt.cpp
  incflo_regrid.cpp

  CFDSim.cpp
  )

target_sources(${amr_wind_exe_name} PRIVATE main.cpp)

#AMR-Wind include directories
target_include_directories(${amr_wind_lib_name} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
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
add_subdirectory(prob)
add_subdirectory(wind_energy)
add_subdirectory(equation_systems)
add_subdirectory(transport_models)
add_subdirectory(turbulence)
add_subdirectory(physics)
add_subdirectory(overset)

if(AMR_WIND_ENABLE_MASA)
  add_subdirectory(mms)
endif()

# Generate AMReX build information file
include(AMReXBuildInfo)
# Create an internal target as amrwind library is an OBJECT library and cannot have
# PRE_BUILD etc.
add_library(amrwind_buildinfo)
generate_buildinfo(amrwind_buildinfo ${CMAKE_SOURCE_DIR})

target_link_libraries_system(${amr_wind_lib_name} PUBLIC AMReX::amrex)
if (AMR_WIND_ENABLE_NETCDF)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_NETCDF)
  target_link_libraries_system(${amr_wind_lib_name} PUBLIC NetCDF::NetCDF)
endif()
target_link_libraries(${amr_wind_lib_name} PRIVATE amrwind_buildinfo)
target_link_libraries(${amr_wind_exe_name} PRIVATE ${amr_wind_lib_name})
target_link_libraries(${aw_api_lib} PUBLIC ${amr_wind_lib_name})

# Set -fpic options
set_target_properties(
  ${amr_wind_lib_name} amrwind_buildinfo
  PROPERTIES POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE amr_wind_header_list ${CMAKE_CURRENT_SOURCE_DIR} "*.H")
add_custom_target(amr_wind_headers SOURCES ${amr_wind_header_list})
