
# target_link_libraries_system
#
# This function is similar to target_link_libraries but allows the includes
# determined from the library to be added as system includes to suppress
# warnings generated from those header files
#
# https://stackoverflow.com/questions/52135983/cmake-target-link-libraries-include-as-system-to-suppress-compiler-warnings/52136398#52136398
#
function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

macro(init_amrex)
  if (${AMR_WIND_USE_INTERNAL_AMREX})
    set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/submods/amrex")
    include(${CMAKE_SOURCE_DIR}/cmake/set_amrex_options.cmake)
    list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
    add_subdirectory(${AMREX_SUBMOD_LOCATION})
    set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/submods/amrex/Tools/Plotfile/fcompare
      CACHE INTERNAL "Path to fcompare executable for regression tests")
  else()
    set(CMAKE_PREFIX_PATH ${AMREX_DIR} ${CMAKE_PREFIX_PATH})
    list(APPEND AMREX_COMPONENTS
      "3D" "PIC" "PARTICLES" "DPARTICLES" "DP" "LSOLVERS")
    if (AMR_WIND_ENABLE_MPI)
      list(APPEND AMREX_COMPONENTS "MPI")
    endif()
    if (AMR_WIND_ENABLE_OPENMP)
      list(APPEND AMREX_COMPONENTS "OMP")
    endif()
    if (AMR_WIND_ENABLE_CUDA)
      list(APPEND AMREX_COMPONENTS "CUDA")
    endif()
    separate_arguments(AMREX_COMPONENTS)
    find_package(AMReX CONFIG REQUIRED
      COMPONENTS ${AMREX_COMPONENTS})
    message(STATUS "Found AMReX = ${AMReX_DIR}")
    set(FCOMPARE_EXE ${AMReX_DIR}/../../../bin/fcompare
      CACHE INTERNAL "Path to fcompare executable for regression tests")
  endif()
endmacro(init_amrex)
