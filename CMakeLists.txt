############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(AMR-Wind CXX C Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for the project
option(AMR_WIND_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(AMR_WIND_ENABLE_SPHINX_API_DOCS "Link Doxygen API docs to Sphinx" OFF)
option(AMR_WIND_ENABLE_ALL_WARNINGS "Show most warnings for most compilers" OFF)
option(AMR_WIND_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)

#Enabling tests overrides the executable options
option(AMR_WIND_ENABLE_TESTS "Enable testing suite" OFF)
option(AMR_WIND_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)

#Options for the executable
option(AMR_WIND_ENABLE_MPI "Enable MPI" OFF)
option(AMR_WIND_ENABLE_OPENMP "Enable OpenMP" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

########################### AMReX #####################################

set(AMREX_SUBMOD_LOCATION "submods/amrex")
include(${CMAKE_SOURCE_DIR}/cmake/set_amrex_options.cmake)
add_subdirectory(${AMREX_SUBMOD_LOCATION})

if(AMR_WIND_ENABLE_FCOMPARE OR AMR_WIND_TEST_WITH_FCOMPARE)
  add_subdirectory(${AMREX_SUBMOD_LOCATION}/Tools/Plotfile)
endif()

########################### AMR-Wind #####################################

find_package(Python REQUIRED)
find_package(Threads REQUIRED) # Needed this for the Travis CI system
if(AMR_WIND_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "AMR-Wind Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# When building, don't use the install RPATH already (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

#Create target names
set(amr_wind_exe_name "amr_wind")
set(amr_wind_unit_test_exe_name "${amr_wind_exe_name}_unit_tests")

#Create main target executable
add_executable(${amr_wind_exe_name} "")

include(${CMAKE_SOURCE_DIR}/cmake/set_compile_flags.cmake)

#Keep our Fortran module files confined to a unique directory for each executable 
set_target_properties(${amr_wind_exe_name} PROPERTIES Fortran_MODULE_DIRECTORY
                     "${CMAKE_BINARY_DIR}/fortran_modules/${amr_wind_exe_name}_fortran_modules")

#Create directory unique to executable to store generated files
set(GENERATED_FILES_DIR ${CMAKE_BINARY_DIR}/generated_files/${amr_wind_exe_name}_generated_files)
file(MAKE_DIRECTORY ${GENERATED_FILES_DIR})

#Generate the AMReX_buildInfo.cpp file with Python
add_custom_target(generate_build_info_${amr_wind_exe_name} ALL
   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}/Tools/C_scripts/makebuildinfo_C.py
   --amrex_home "${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}"
   --COMP ${CMAKE_CXX_COMPILER_ID} --COMP_VERSION ${CMAKE_CXX_COMPILER_VERSION}
   --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_Fortran_COMPILER_VERSION}
   --GIT "${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}"
   WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
   COMMENT "Generating AMReX_buildInfo.cpp"
)                  
  
#Set the dependencies on targets so the generated source code files are there before we try to build the executable 
add_dependencies(${amr_wind_exe_name} generate_build_info_${amr_wind_exe_name})

#Build amr-wind and link to amrex library
add_subdirectory(src)

if(AMR_WIND_ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif()

if(AMR_WIND_ENABLE_DOCUMENTATION)
   add_subdirectory(docs)
endif()

#Define what we want to be installed during a make install 
install(TARGETS ${amr_wind_exe_name}
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
