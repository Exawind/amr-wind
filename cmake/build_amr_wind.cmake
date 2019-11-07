function(build_amr_wind amr_wind_exe_name amr_wind_exe_options_file)

  unset(AMR_WIND_EXTRA_SOURCES)
  unset(AMR_WIND_DIM)
  unset(AMR_WIND_ENABLE_EB)
  unset(AMR_WIND_ENABLE_MASA)

  get_filename_component(exe_directory ${amr_wind_exe_options_file} DIRECTORY)
  include(${amr_wind_exe_options_file})

  #message("-- AMR_WIND_DIM = ${AMR_WIND_DIM}D")

  if(AMR_WIND_ENABLE_EB)
    set(EB "eb")
  else()
    unset(EB)
  endif()

  #Expose functions we want to be able to call
  include(${CMAKE_SOURCE_DIR}/cmake/add_source_function.cmake)
  include(${CMAKE_SOURCE_DIR}/cmake/amr_wind_sources.cmake)

  #Clear source file list from any previous executables
  set_property(GLOBAL PROPERTY GlobalSourceList "") 

  #Aggregate amrex source files
  get_amr_wind_sources(${amr_wind_exe_name})
  
  #Put source list from global property into local list
  get_property(AMR_WIND_SOURCES GLOBAL PROPERTY GlobalSourceList)

  #Create the full path to the extra case-specific source files
  #Each AMR_WIND_EXTRA_SOURCE must be an explicit path to each source file at the moment
  foreach(AMR_WIND_EXTRA_SOURCE ${AMR_WIND_EXTRA_SOURCES})
    list(APPEND MY_EXTRA_SOURCES ${AMR_WIND_EXTRA_SOURCE})
  endforeach()

  #Create an executable based on all the source files we aggregated
  add_executable(${amr_wind_exe_name} ${AMR_WIND_SOURCES} ${MY_EXTRA_SOURCES})
  target_link_libraries(${amr_wind_exe_name} PRIVATE amrex${AMR_WIND_DIM}d${EB})

  #AMReX definitions
  target_compile_definitions(${amr_wind_exe_name} PRIVATE BL_SPACEDIM=${AMR_WIND_DIM})
  target_compile_definitions(${amr_wind_exe_name} PRIVATE BL_FORT_USE_UNDERSCORE)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_SPACEDIM=${AMR_WIND_DIM})
  target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_FORT_USE_UNDERSCORE)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
  target_compile_definitions(${amr_wind_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)

  # CMake BUILD_TYPE should already define this
  #if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  #  target_compile_definitions(${amr_wind_exe_name} PRIVATE NDEBUG)
  #endif()

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
  target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/src/diffusion)
  target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/src/rheology)
  target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/src/boundary_conditions)
  target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/src/setup)
  target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/src/derive)
  if(AMR_WIND_ENABLE_EB)
    target_include_directories(${amr_wind_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/src/embedded_boundaries)
  endif()
  
  #Link our executable to the MPI libraries, etc
  if(AMR_WIND_ENABLE_MPI)
    target_link_libraries(${amr_wind_exe_name} PRIVATE MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran)
    target_compile_definitions(${amr_wind_exe_name} PRIVATE BL_USE_MPI)
    target_compile_definitions(${amr_wind_exe_name} PRIVATE AMREX_USE_MPI)
  endif()

  #Link our executable to the MASA libraries, etc
  if(AMR_WIND_ENABLE_MASA)
    target_link_libraries(${amr_wind_exe_name} PRIVATE ${MASA_LIBRARY} ${MASA_FORTRAN_LIBRARY})
    target_compile_definitions(${amr_wind_exe_name} PRIVATE USE_MASA)
    target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${MASA_INCLUDE_DIRS})
    target_include_directories(${amr_wind_exe_name} SYSTEM PRIVATE ${MASA_MOD_DIRS})
  endif()

  #if(AMR_WIND_ENABLE_OPENMP)
  #  string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
  #  string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
  #  string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
  #endif()

  #Keep our Fortran module files confined to a unique directory for each executable 
  set_target_properties(${amr_wind_exe_name} PROPERTIES Fortran_MODULE_DIRECTORY
                       "${CMAKE_BINARY_DIR}/fortran_modules/${amr_wind_exe_name}_fortran_modules")

  #Create directory unique to executable to store generated files
  set(GENERATED_FILES_DIR ${CMAKE_BINARY_DIR}/generated_files/${amr_wind_exe_name}_generated_files)
  file(MAKE_DIRECTORY ${GENERATED_FILES_DIR})

  #set(PARAMETER_DIRS "")
  #string(APPEND PARAMETER_DIRS " ${CMAKE_SOURCE_DIR}/submods/PelePhysics/Eos/${AMR_WIND_EOS_MODEL}/_parameters")

  if(PYTHON_FOUND)
     #Generate the extern.f90 file with Python
     #add_custom_target(generate_extern_${amr_wind_exe_name} ALL
     #   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/submods/amrex/Tools/F_scripts/write_probin.py
     #   -t ${CMAKE_SOURCE_DIR}/src/extern_probin.template
     #   -o extern.f90 -n extern --pa "${PARAMETER_DIRS}"
     #   WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/extern.f90
     #   COMMENT "Generating extern.f90"
     #)

     #Generate the AMReX_buildInfo.cpp file with Python
     add_custom_target(generate_build_info_${amr_wind_exe_name} ALL
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/submods/amrex/Tools/C_scripts/makebuildinfo_C.py
        --amrex_home "${CMAKE_SOURCE_DIR}/submods/amrex"                        
        --COMP ${CMAKE_C_COMPILER_ID} --COMP_VERSION ${CMAKE_C_COMPILER_VERSION}
        --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_Fortran_COMPILER_VERSION}
        --GIT "${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/submods/amrex"
        WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
        COMMENT "Generating AMReX_buildInfo.cpp"
     )                  
  endif()
  
  #Set the dependencies on targets so the generated source code files are there before we try to build the executable 
  #add_dependencies(${amr_wind_exe_name} generate_extern_${amr_wind_exe_name} generate_build_info_${amr_wind_exe_name})
  add_dependencies(${amr_wind_exe_name} generate_build_info_${amr_wind_exe_name})
 
  #Define what we want to be installed during a make install 
  install(TARGETS ${amr_wind_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
  
endfunction(build_amr_wind amr_wind_exe_name amr_wind_exe_options_file)
