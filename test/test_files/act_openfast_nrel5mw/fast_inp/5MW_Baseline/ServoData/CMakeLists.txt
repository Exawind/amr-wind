# Need to set -w (no warnings) for GCC < 7 so it doesn't fail when trying to output
# warnings. This problem was fixed in GCC 7.
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-w ${CMAKE_Fortran_FLAGS}")
endif()
add_library(DISCON SHARED DISCON.F90)
