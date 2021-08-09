############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(AMR-Wind CXX C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(amr-wind-utils)

########################## OPTIONS #####################################

#General options for the project
option(AMR_WIND_ENABLE_ALL_WARNINGS "Show most warnings for most compilers" ON)
option(AMR_WIND_ENABLE_WERROR "Treat compiler warnings as errors" OFF)
option(AMR_WIND_ENABLE_CLANG_TIDY "Compile with clang-tidy static analysis" OFF)
option(AMR_WIND_ENABLE_CPPCHECK "Enable cppcheck static analysis target" OFF)
option(AMR_WIND_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)

#Enabling tests overrides the executable options
option(AMR_WIND_ENABLE_UNIT_TESTS "Enable unit testing" ON)
option(AMR_WIND_ENABLE_TESTS "Enable testing suite" OFF)
option(AMR_WIND_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)
option(AMR_WIND_SAVE_GOLDS "Provide a directory in which to save golds during testing" OFF)

#Options for the executable
option(AMR_WIND_ENABLE_MPI "Enable MPI" OFF)
option(AMR_WIND_ENABLE_OPENMP "Enable OpenMP" OFF)
option(AMR_WIND_ENABLE_CUDA "Enable CUDA" OFF)
option(AMR_WIND_ENABLE_HIP "Enable HIP/ROCm" OFF)
option(AMR_WIND_ENABLE_DPCPP "Enable Intel OneAPI DPC++" OFF)
option(AMR_WIND_ENABLE_TINY_PROFILE "Enable AMReX TinyProfile support" OFF)

# Third party libraries
option(AMR_WIND_USE_INTERNAL_AMREX "Add AMReX as subproject" ON)
option(AMR_WIND_ENABLE_NETCDF "Enable NetCDF library" OFF)
option(AMR_WIND_ENABLE_MASA "Enable MASA library" OFF)
option(AMR_WIND_ENABLE_HYPRE "Enable HYPRE integration" OFF)
option(AMR_WIND_ENABLE_OPENFAST "Enable OpenFAST integration" OFF)
option(AMR_WIND_ENABLE_ASCENT "Enable Ascent visualization library" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(AMR_WIND_ENABLE_CUDA)
  enable_language(CUDA)
  # Fix issues with GPU builds
  enable_language(Fortran)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
endif()

if (AMR_WIND_ENABLE_HIP)
  find_package(HIP REQUIRED)
  if (NOT DEFINED AMD_ARCH)
    # Set default AMD architectures (based on Tulip for now)
    set(AMD_ARCH "gfx906,gfx908")
  endif()
endif()

if(AMR_WIND_TEST_WITH_FCOMPARE)
  set(AMR_WIND_ENABLE_FCOMPARE ON)
endif()

if (AMR_WIND_ENABLE_OPENFAST)
  enable_language(Fortran)
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

include(set_rpath)

#Create target names
set(amr_wind_lib_name "amrwind_obj")
set(amr_wind_exe_name "amr_wind")
set(amr_wind_unit_test_exe_name "${amr_wind_exe_name}_unit_tests")
set(aw_api_lib "amrwind_api")

#Create main target executable
add_library(${amr_wind_lib_name} OBJECT)
add_library(${aw_api_lib} SHARED)
add_executable(${amr_wind_exe_name})

init_code_checks()
if(CLANG_TIDY_EXE)
  set_target_properties(${amr_wind_lib_name} ${aw_api_lib} ${amr_wind_exe_name}
                        PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
endif()

include(set_compile_flags)

if (AMR_WIND_ENABLE_NETCDF)
  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
  find_package(NetCDF QUIET REQUIRED)
  if(NetCDF_FOUND)
    message(STATUS "Found NetCDF = ${NETCDF_DIR}")
  endif()
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_NETCDF)
  target_link_libraries_system(${amr_wind_lib_name} PUBLIC NetCDF::NetCDF)
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

if (AMR_WIND_ENABLE_OPENFAST)
  set(CMAKE_PREFIX_PATH ${OPENFAST_DIR} ${CMAKE_PREFIX_PATH})
  find_package(OpenFAST QUIET REQUIRED)
  if(OpenFAST_FOUND)
    message(STATUS "Found OpenFAST = ${OPENFAST_DIR}")
  endif()
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_OPENFAST)
  target_include_directories(${amr_wind_lib_name} PUBLIC ${OpenFAST_INCLUDE_DIRS})
  target_link_libraries(${amr_wind_lib_name} PUBLIC ${OpenFAST_LIBRARIES})
endif()

if(AMR_WIND_ENABLE_ASCENT)
  set(CMAKE_PREFIX_PATH ${ASCENT_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Ascent QUIET REQUIRED)
  if (Ascent_FOUND)
    message(STATUS "Found Ascent = ${ASCENT_INCLUDE_DIRS}")
  endif()
  target_link_libraries(${amr_wind_lib_name} PUBLIC ascent::ascent)
  target_compile_definitions(${amr_wind_lib_name} PRIVATE AMR_WIND_USE_ASCENT)
endif()

#Build amr-wind and link to amrex library
generate_version_info()
add_subdirectory(amr-wind)

if (AMR_WIND_ENABLE_CUDA)
  set(awtargets "${amr_wind_lib_name};${amr_wind_exe_name}")
  foreach(tgt IN LISTS awtargets)
    set_cuda_build_properties(${tgt})
    # get_target_property(AMR_WIND_SOURCES ${tgt} SOURCES)
    # list(FILTER AMR_WIND_SOURCES INCLUDE REGEX "\\.cpp")
    # set_source_files_properties(${AMR_WIND_SOURCES} PROPERTIES LANGUAGE CUDA)
  endforeach()
endif()

if (AMR_WIND_ENABLE_UNIT_TESTS OR AMR_WIND_ENABLE_TESTS)
  add_executable(${amr_wind_unit_test_exe_name})
  if(CLANG_TIDY_EXE)
    set_target_properties(${amr_wind_unit_test_exe_name}
                          PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
  endif()
  add_subdirectory("submods/googletest")
  add_subdirectory("unit_tests")
  set_cuda_build_properties(${amr_wind_unit_test_exe_name})
  # if (AMR_WIND_ENABLE_CUDA)
  #   get_target_property(UTEST_SOURCES ${amr_wind_unit_test_exe_name} SOURCES)
  #   set_source_files_properties(${UTEST_SOURCES} PROPERTIES LANGUAGE CUDA)
  #   set_target_properties(${amr_wind_unit_test_exe_name} PROPERTIES
  #     CUDA_SEPARABLE_COMPILATION ON)
  # endif()
endif()

add_subdirectory(tools)

if(AMR_WIND_ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif()

#Define what we want to be installed during a make install
install(
  TARGETS ${amr_wind_exe_name} ${amr_wind_unit_test_exe_name}
  ${aw_api_lib} ${amr_wind_lib_name} buildInfo${amr_wind_lib_name}
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
