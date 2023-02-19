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
option(AMR_WIND_ENABLE_CLANG_TIDY "Compile with clang-tidy static analysis" OFF)
option(AMR_WIND_ENABLE_CPPCHECK "Enable cppcheck static analysis target" OFF)
option(AMR_WIND_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)

#Enabling tests overrides the executable options
option(AMR_WIND_ENABLE_UNIT_TESTS "Enable unit testing" ON)
option(AMR_WIND_ENABLE_TESTS "Enable testing suite" OFF)
option(AMR_WIND_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)
option(AMR_WIND_SAVE_GOLDS "Provide a directory in which to save golds during testing" OFF)
option(AMR_WIND_ENABLE_FPE_TRAP_FOR_TESTS "Enable FPE trapping in tests" ON)

#Options for the executable
option(AMR_WIND_ENABLE_MPI "Enable MPI" OFF)
option(AMR_WIND_ENABLE_OPENMP "Enable OpenMP" OFF)
option(AMR_WIND_ENABLE_CUDA "Enable CUDA" OFF)
option(AMR_WIND_ENABLE_ROCM "Enable ROCm/HIP" OFF)
option(AMR_WIND_ENABLE_DPCPP "Enable Intel OneAPI DPC++" OFF)
option(AMR_WIND_ENABLE_TINY_PROFILE "Enable AMReX TinyProfile support" OFF)
set(AMR_WIND_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

# Third party libraries
option(AMR_WIND_USE_INTERNAL_AMREX "Add AMReX as subproject" ON)
option(AMR_WIND_USE_INTERNAL_AMREX_HYDRO "Add AMReX-Hydro as subproject" ON)
option(AMR_WIND_ENABLE_HDF5 "Enable HDF5 library" OFF)
option(AMR_WIND_ENABLE_HDF5_ZFP "Enable ZFP compression in HDF5 library" OFF)
option(AMR_WIND_ENABLE_NETCDF "Enable NetCDF library" OFF)
option(AMR_WIND_ENABLE_MASA "Enable MASA library" OFF)
option(AMR_WIND_ENABLE_HELICS "Enable HELICS library" OFF)
option(AMR_WIND_ENABLE_HYPRE "Enable HYPRE integration" OFF)
option(AMR_WIND_ENABLE_OPENFAST "Enable OpenFAST integration" OFF)
option(AMR_WIND_ENABLE_ASCENT "Enable Ascent visualization library" OFF)
option(AMR_WIND_ENABLE_UMPIRE "Enable Umpire GPU memory pools" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(AMR_WIND_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
  # Fix issues with GPU builds
  if(AMR_WIND_ENABLE_OPENFAST)
    enable_language(Fortran)
  endif()
endif()

if(AMR_WIND_ENABLE_ROCM)
  find_package(HIP REQUIRED)
  if(NOT DEFINED AMD_ARCH)
    # Set default AMD architectures (based on Frontier)
    set(AMD_ARCH "gfx90a")
  endif()
endif()

if (AMR_WIND_ENABLE_UMPIRE)
  set(CMAKE_PREFIX_PATH ${UMPIRE_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Umpire REQUIRED)
  include_directories(SYSTEM ${UMPIRE_INCLUDE_DIRS})
endif()

if(AMR_WIND_TEST_WITH_FCOMPARE)
  set(AMR_WIND_ENABLE_FCOMPARE ON)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(AMR_WIND_ENABLE_FPE_TRAP_FOR_TESTS OFF)
  message(WARNING "Disabling FPE trapping for tests when using AppleClang.")
endif()

########################### AMReX #####################################
init_amrex()
init_amrex_hydro()

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
add_library(${aw_api_lib})
add_executable(${amr_wind_exe_name})

init_code_checks()
if(CLANG_TIDY_EXE)
  set_target_properties(${amr_wind_lib_name} ${aw_api_lib} ${amr_wind_exe_name}
                        PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
endif()

include(set_compile_flags)

if(AMR_WIND_ENABLE_NETCDF)
  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
  find_package(NetCDF REQUIRED)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_NETCDF)
  target_link_libraries_system(${amr_wind_lib_name} PUBLIC NetCDF::NetCDF)
endif()

if(AMR_WIND_ENABLE_HDF5)
  set(CMAKE_PREFIX_PATH ${HDF5_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_HDF5)
endif()

if(AMR_WIND_ENABLE_HDF5_ZFP)
  if(NOT AMR_WIND_ENABLE_HDF5)
    message(FATAL_ERROR "HDF5 must be enabled to use ZFP compression in HDF5 plot files")
  endif()
  set(CMAKE_PREFIX_PATH ${H5Z_ZFP_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_HDF5_ZFP)
endif()

if(AMR_WIND_ENABLE_MASA)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA REQUIRED)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_MASA)
  target_link_libraries(${amr_wind_lib_name} PUBLIC MASA::MASA)
endif()

if(AMR_WIND_ENABLE_OPENFAST)
  set(CMAKE_PREFIX_PATH ${OPENFAST_DIR} ${CMAKE_PREFIX_PATH})
  find_package(OpenFAST REQUIRED)
  include(CheckIncludeFileCXX)
  set(CMAKE_REQUIRED_INCLUDES "${OpenFAST_INCLUDE_DIRS}")
  check_include_file_cxx("SCDataEx_Types.h" USE_SCDX)
  if(USE_SCDX)
    target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_FAST_USE_SCDX)
  endif()
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_OPENFAST)
  target_include_directories(${amr_wind_lib_name} PUBLIC ${OpenFAST_INCLUDE_DIRS})
  target_link_libraries(${amr_wind_lib_name} PUBLIC ${OpenFAST_LIBRARIES})
endif()

if(AMR_WIND_ENABLE_ASCENT)
  set(CMAKE_PREFIX_PATH ${ASCENT_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Ascent REQUIRED)
  if(AMR_WIND_ENABLE_MPI)
    target_link_libraries(${amr_wind_lib_name} PUBLIC ascent::ascent_mpi)
  else()
    target_link_libraries(${amr_wind_lib_name} PUBLIC ascent::ascent)
  endif()
  target_compile_definitions(${amr_wind_lib_name} PRIVATE AMR_WIND_USE_ASCENT)
endif()

# Link with HELICS module 
if(AMR_WIND_ENABLE_HELICS)
  set(CMAKE_PREFIX_PATH ${HELICS_DIR} ${CMAKE_PREFIX_PATH})
  find_package(HELICS 3 REQUIRED)
  target_include_directories(${amr_wind_lib_name} PUBLIC ${HELICS_INSTALL_DIR}/include)
  target_link_libraries(${amr_wind_lib_name} PUBLIC HELICS::helicsCpp98)
  target_compile_definitions(${amr_wind_lib_name} PRIVATE AMR_WIND_USE_HELICS)
endif()

#Build amr-wind and link to amrex library
generate_version_info()
add_subdirectory(amr-wind)

if(AMR_WIND_ENABLE_CUDA)
  set(awtargets "${amr_wind_lib_name};${amr_wind_exe_name}")
  foreach(tgt IN LISTS awtargets)
    set_cuda_build_properties(${tgt})
    # get_target_property(AMR_WIND_SOURCES ${tgt} SOURCES)
    # list(FILTER AMR_WIND_SOURCES INCLUDE REGEX "\\.cpp")
    # set_source_files_properties(${AMR_WIND_SOURCES} PROPERTIES LANGUAGE CUDA)
  endforeach()
endif()

if(AMR_WIND_ENABLE_UNIT_TESTS OR AMR_WIND_ENABLE_TESTS)
  add_executable(${amr_wind_unit_test_exe_name})
  if(CLANG_TIDY_EXE)
    set_target_properties(${amr_wind_unit_test_exe_name}
                          PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
  endif()
  add_subdirectory("submods/googletest")
  add_subdirectory("unit_tests")
  set_cuda_build_properties(${amr_wind_unit_test_exe_name})
  # if(AMR_WIND_ENABLE_CUDA)
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
  TARGETS ${amr_wind_exe_name}
  ${aw_api_lib} ${amr_wind_lib_name} buildInfo${amr_wind_lib_name}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

if(AMR_WIND_ENABLE_UNIT_TESTS OR AMR_WIND_ENABLE_TESTS)
  install(
    TARGETS ${amr_wind_unit_test_exe_name}
    RUNTIME DESTINATION bin)
endif()

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
