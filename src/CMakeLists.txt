target_sources(${amr_wind_exe_name}
   PRIVATE
      #C++
      main.cpp
      advance.cpp
      incflo.cpp
      incflo_compute_dt.cpp
      #Headers
      incflo_bclist.H
      param_mod_F.H
      incflo_proj_F.H
      #Fortran
      param_mod.f90
   )

# Set library suffixes for EB enabled
if(AMR_WIND_ENABLE_EB)
  set(EB "eb")
else()
  unset(EB)
endif()

target_link_libraries(${amr_wind_exe_name} PRIVATE amrex${AMR_WIND_DIM}d${EB})

#AMReX definitions
target_compile_definitions(${amr_wind_exe_name} PRIVATE BL_SPACEDIM=${AMR_WIND_DIM})
target_compile_definitions(${amr_wind_exe_name} PRIVATE BL_FORT_USE_UNDERSCORE)
target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_SPACEDIM=${AMR_WIND_DIM})
target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_FORT_USE_UNDERSCORE)
target_compile_definitions(${amr_wind_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
target_compile_definitions(${amr_wind_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_compile_definitions(${amr_wind_exe_name} PRIVATE BL_Darwin)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_Darwin)
endif()

#AMR-Wind definitions
if(AMR_WIND_ENABLE_EB)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE AMR_WIND_USE_EB)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_USE_EB)
endif()

#AMReX include directories
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/Base)
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/AmrCore)
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/Boundary)
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/LinearSolvers/MLMG)
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/LinearSolvers/Projections)
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Tools/C_scripts)
target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_BINARY_DIR}/fortran_modules/amrex${AMR_WIND_DIM}d${EB}_fortran_modules)
if(AMR_WIND_ENABLE_EB)
  target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/EB)
endif()

#AMR-Wind include directories
target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  
#Link our executable to the MPI libraries, etc
if(AMR_WIND_ENABLE_MPI)
  target_link_libraries(${amr_wind_exe_name} PRIVATE MPI::MPI_CXX MPI::MPI_Fortran)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE BL_USE_MPI)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_USE_MPI)
endif()

#Add generated source files
target_sources(${amr_wind_exe_name}
   PRIVATE
      ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
)
 
add_subdirectory(boundary_conditions)
add_subdirectory(convection)
add_subdirectory(derive)
add_subdirectory(diffusion)
add_subdirectory(projection)
add_subdirectory(rheology)
add_subdirectory(setup)
add_subdirectory(utilities)
if(AMR_WIND_ENABLE_EB)
   add_subdirectory(embedded_boundaries)
endif()
