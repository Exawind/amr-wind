add_executable(${amr_wind_unit_test_exe_name} "")

target_sources(${amr_wind_unit_test_exe_name}
   PRIVATE
      #C++
      unit_tests_3d.C
      unit_test_3d_1.C
   )

# Set library suffixes for EB enabled
if(AMR_WIND_ENABLE_EB)
  set(EB "eb")
else()
  unset(EB)
endif()

target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE amrex${AMR_WIND_DIM}d${EB})

#AMReX definitions
target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE BL_SPACEDIM=${AMR_WIND_DIM})
target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE BL_FORT_USE_UNDERSCORE)
target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE AMREX_SPACEDIM=${AMR_WIND_DIM})
target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE AMREX_FORT_USE_UNDERSCORE)
target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE BL_Darwin)
  target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE AMREX_Darwin)
endif()

#AMR-Wind definitions
if(AMR_WIND_ENABLE_EB)
  target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE AMR_WIND_USE_EB)
  target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE AMREX_USE_EB)
endif()

#AMReX include directories
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/Base)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/AmrCore)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/Boundary)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/LinearSolvers/MLMG)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/LinearSolvers/Projections)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Tools/C_scripts)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_BINARY_DIR}/fortran_modules/amrex${AMR_WIND_DIM}d${EB}_fortran_modules)
if(AMR_WIND_ENABLE_EB)
  target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/amrex/Src/EB)
endif()

#AMR-Wind include directories
target_include_directories(${amr_wind_unit_test_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  
#Link our executable to the MPI libraries, etc
if(AMR_WIND_ENABLE_MPI)
  target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE MPI::MPI_CXX MPI::MPI_Fortran)
  target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE BL_USE_MPI)
  target_compile_definitions(${amr_wind_unit_test_exe_name} PRIVATE AMREX_USE_MPI)
endif()

#Link our unit test executable with GoogleTest
#target_link_libraries(${utest_exe_name} PRIVATE GTest::GTest GTest::Main)
target_link_libraries(${amr_wind_unit_test_exe_name} PRIVATE gtest gtest_main)
target_include_directories(${amr_wind_unit_test_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/submods/googletest/googletest/include)

#Keep our Fortran module files confined to a unique directory for each executable 
set_target_properties(${amr_wind_unit_test_exe_name} PROPERTIES Fortran_MODULE_DIRECTORY
                     "${CMAKE_BINARY_DIR}/fortran_modules/${amr_wind_unit_test_exe_name}_fortran_modules")

#Create directory unique to executable to store generated files
set(GENERATED_FILES_DIR ${CMAKE_BINARY_DIR}/generated_files/${amr_wind_unit_test_exe_name}_generated_files)
file(MAKE_DIRECTORY ${GENERATED_FILES_DIR})

#Generate the AMReX_buildInfo.cpp file with Python
add_custom_target(generate_build_info_${amr_wind_unit_test_exe_name} ALL
   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/submods/amrex/Tools/C_scripts/makebuildinfo_C.py
   --amrex_home "${CMAKE_SOURCE_DIR}/submods/amrex"                        
   --COMP ${CMAKE_CXX_COMPILER_ID} --COMP_VERSION ${CMAKE_CXX_COMPILER_VERSION}
   --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_Fortran_COMPILER_VERSION}
   --GIT "${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/submods/amrex"
   WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
   COMMENT "Generating AMReX_buildInfo.cpp"
)                  
  
#Set the dependencies on targets so the generated source code files are there before we try to build the executable 
add_dependencies(${amr_wind_unit_test_exe_name} generate_build_info_${amr_wind_unit_test_exe_name})

#Add generated source files
target_sources(${amr_wind_unit_test_exe_name}
   PRIVATE
      ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
)
 
#Define what we want to be installed during a make install 
install(TARGETS ${amr_wind_unit_test_exe_name}
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
