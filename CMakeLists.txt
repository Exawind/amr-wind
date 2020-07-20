############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(AMR-Wind CXX C)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(${CMAKE_SOURCE_DIR}/cmake/amr-wind-utils.cmake)

########################## OPTIONS #####################################

#General options for the project
option(AMR_WIND_ENABLE_FORTRAN "Enable Fortran interfaces" OFF)
option(AMR_WIND_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(AMR_WIND_ENABLE_SPHINX_API_DOCS "Link Doxygen API docs to Sphinx" OFF)
option(AMR_WIND_ENABLE_ALL_WARNINGS "Show most warnings for most compilers" OFF)
option(AMR_WIND_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)

#Enabling tests overrides the executable options
option(AMR_WIND_ENABLE_UNIT_TESTS "Enable unit testing" ON)
option(AMR_WIND_ENABLE_TESTS "Enable testing suite" OFF)
option(AMR_WIND_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)

#Options for the executable
option(AMR_WIND_ENABLE_MPI "Enable MPI" OFF)
option(AMR_WIND_ENABLE_OPENMP "Enable OpenMP" OFF)
option(AMR_WIND_ENABLE_CUDA "Enable CUDA" OFF)
option(AMR_WIND_ENABLE_DPCPP "Enable Intel OneAPI DPC++" OFF)

# Third party libraries
option(AMR_WIND_USE_INTERNAL_AMREX "Add AMReX as subproject" ON)
option(AMR_WIND_ENABLE_NETCDF "Enable NetCDF library" OFF)
option(AMR_WIND_ENABLE_MASA "Enable MASA library" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(AMR_WIND_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
endif()

if(AMR_WIND_TEST_WITH_FCOMPARE)
  set(AMR_WIND_ENABLE_FCOMPARE ON)
endif()

########################### AMReX #####################################
init_amrex()

########################### AMR-Wind #####################################

if(AMR_WIND_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "AMR-Wind Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

include(${CMAKE_SOURCE_DIR}/cmake/set_rpath.cmake)

#Create target names
set(amr_wind_lib_name "amrwind_obj")
set(amr_wind_exe_name "amr_wind")
set(amr_wind_unit_test_exe_name "${amr_wind_exe_name}_unit_tests")
set(aw_api_lib "amrwind_api")

#Create main target executable
add_library(${amr_wind_lib_name} OBJECT)
add_library(${aw_api_lib} SHARED)
add_executable(${amr_wind_exe_name})

include(${CMAKE_SOURCE_DIR}/cmake/set_compile_flags.cmake)

if (AMR_WIND_ENABLE_FORTRAN)
  #Keep our Fortran module files confined to a directory
  set_target_properties(
    ${amr_wind_lib_name} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran_modules/")
endif()

if (AMR_WIND_ENABLE_NETCDF)
  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
  find_package(NetCDF REQUIRED)
endif()

if(AMR_WIND_ENABLE_MASA)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA QUIET REQUIRED)
  if(MASA_FOUND)
    message(STATUS "Found MASA = ${MASA_DIR}")
  endif()
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_MASA)
  target_link_libraries(${amr_wind_lib_name} PUBLIC MASA::MASA)
endif()

#Build amr-wind and link to amrex library
add_subdirectory(amr-wind)

if (AMR_WIND_ENABLE_CUDA)
  set(awtargets "${amr_wind_lib_name};${amr_wind_exe_name}")
  foreach(tgt IN LISTS awtargets)
    get_target_property(AMR_WIND_SOURCES ${tgt} SOURCES)
    list(FILTER AMR_WIND_SOURCES INCLUDE REGEX "\\.cpp")
    set_source_files_properties(${AMR_WIND_SOURCES} PROPERTIES LANGUAGE CUDA)
  endforeach()
endif()

if (AMR_WIND_ENABLE_UNIT_TESTS OR AMR_WIND_ENABLE_TESTS)
  add_executable(${amr_wind_unit_test_exe_name})
  add_subdirectory("submods/googletest")
  add_subdirectory("unit_tests")
  if (AMR_WIND_ENABLE_CUDA)
    get_target_property(UTEST_SOURCES ${amr_wind_unit_test_exe_name} SOURCES)
    set_source_files_properties(${UTEST_SOURCES} PROPERTIES LANGUAGE CUDA)
  endif()
endif()

if(AMR_WIND_ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif()

if(AMR_WIND_ENABLE_DOCUMENTATION)
   add_subdirectory(docs)
endif()

#Define what we want to be installed during a make install
install(
  TARGETS ${amr_wind_exe_name} ${amr_wind_unit_test_exe_name}
  ${aw_api_lib} ${amr_wind_lib_name} amrwind_buildinfo
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/amr-wind
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.H")

install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION lib/cmake/${PROJECT_NAME})

configure_package_config_file(
  cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION lib/cmake/${PROJECT_NAME})
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})
install(FILES
  ${PROJECT_SOURCE_DIR}/cmake/FindMASA.cmake
  ${PROJECT_SOURCE_DIR}/cmake/FindNetCDF.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Modules)
