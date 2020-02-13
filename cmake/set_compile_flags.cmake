# Logic for handling warnings
list(APPEND AMR_WIND_CXX_FLAGS "-Wno-pass-failed") # Ignore loop not vectorized warnings
if(AMR_WIND_ENABLE_ALL_WARNINGS)
  # GCC, Clang, and Intel seem to accept these
  list(APPEND AMR_WIND_CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
  if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # ifort doesn't like -Wall
    list(APPEND AMR_WIND_Fortran_FLAGS "-Wall")
  else()
    # Intel always reports some diagnostics we don't necessarily care about
    list(APPEND AMR_WIND_CXX_FLAGS "-diag-disable:11074,11076")
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
    # Avoid notes about -faligned-new with GCC > 7
    list(APPEND AMR_WIND_CXX_FLAGS "-faligned-new")
  endif()
endif()

# Add our extra flags according to language
separate_arguments(AMR_WIND_CXX_FLAGS)
separate_arguments(AMR_WIND_Fortran_FLAGS)
target_compile_options(${amr_wind_exe_name} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${AMR_WIND_CXX_FLAGS}>)
target_compile_options(${amr_wind_exe_name} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${AMR_WIND_Fortran_FLAGS}>)
