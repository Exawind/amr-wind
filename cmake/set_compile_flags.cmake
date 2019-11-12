function(set_compile_flags)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)

  # Add any extra flags based on compiler and/or OS
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
     "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
     "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    # Always add for this compiler
    list(APPEND AMR_WIND_CXX_FLAGS "")
    list(APPEND AMR_WIND_C_FLAGS "")
    list(APPEND AMR_WIND_Fortran_FLAGS "-ffree-line-length-none"
                                    "-ffixed-line-length-none"
                                    "-fno-range-check"
                                    "-fno-second-underscore")
    if(NOT "${BUILD_TYPE}" STREQUAL "DEBUG")
      # Add extra optimization flags
      list(APPEND AMR_WIND_CXX_FLAGS     "")
      list(APPEND AMR_WIND_C_FLAGS       "")
      list(APPEND AMR_WIND_Fortran_FLAGS "")
    else()
      # Add extra debug flags
      list(APPEND AMR_WIND_CXX_FLAGS     "")
      list(APPEND AMR_WIND_C_FLAGS       "")
      list(APPEND AMR_WIND_Fortran_FLAGS "-fcheck=bounds"
                                      "-fbacktrace"
                                      "-ffpe-trap=invalid,zero"
                                      "-finit-real=snan"
                                      "-finit-integer=2147483647"
                                      "-ftrapv")
    endif()
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # Always add for this compiler
    list(APPEND AMR_WIND_CXX_FLAGS     "-restrict")
    list(APPEND AMR_WIND_C_FLAGS       "-restrict")
    list(APPEND AMR_WIND_Fortran_FLAGS "")
    if(NOT "${BUILD_TYPE}" STREQUAL "DEBUG")
      # Add extra optimization flags
      list(APPEND AMR_WIND_CXX_FLAGS     "-ip"
                                      "-qopt-report=5"
                                      "-qopt-report-phase=vec"
                                      "-diag-disable:10397") # Hundreds of remarks about .optrpt generated in the "output" location
      list(APPEND AMR_WIND_C_FLAGS       "-ip"
                                      "-qopt-report=5"
                                      "-qopt-report-phase=vec")
      list(APPEND AMR_WIND_Fortran_FLAGS "-diag-disable:8291") # Remark about high precision in an stdout write statement that a mainframe can't handle
    else()
      # Add extra debug flags
      list(APPEND AMR_WIND_CXX_FLAGS     "-traceback"
                                      "-Wcheck")
      list(APPEND AMR_WIND_C_FLAGS       "-traceback"
                                      "-Wcheck")
      list(APPEND AMR_WIND_Fortran_FLAGS "-traceback"
                                      "-check bounds,uninit,pointers")
    endif()
  endif()
  
  # Logic for managing warnings
  if(ENABLE_ALL_WARNINGS)
    # GCC, Clang, and Intel seem to accept these
    list(APPEND AMR_WIND_CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
    list(APPEND AMR_WIND_C_FLAGS "-Wall" "-Wextra" "-pedantic")
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      # ifort doesn't like -Wall
      list(APPEND AMR_WIND_Fortran_FLAGS "-Wall")
    else()
      # Intel reports some diagnostics we don't necessarily care about
      list(APPEND AMR_WIND_CXX_FLAGS "-diag-disable:11074,11076")
    endif()
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
      # Avoid notes about -faligned-new with GCC > 7
      list(APPEND AMR_WIND_CXX_FLAGS "-faligned-new")
    endif()
  else()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
       "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
      # Avoid warning about explicit vectorization not happening in AMReX
      list(APPEND AMR_WIND_CXX_FLAGS "-Wno-pass-failed")
      list(APPEND AMR_WIND_C_FLAGS "-Wno-pass-failed")
    endif()
  endif()

  # Since we created the flags as lists, we replace the semicolons when it gets converted to a string 
  string(REPLACE ";" " " AMR_WIND_CXX_FLAGS "${AMR_WIND_CXX_FLAGS}")
  string(REPLACE ";" " " AMR_WIND_C_FLAGS "${AMR_WIND_C_FLAGS}")
  string(REPLACE ";" " " AMR_WIND_Fortran_FLAGS "${AMR_WIND_Fortran_FLAGS}")

  # Export regular cmake flags to the parent scope
  set(AMR_WIND_CXX_FLAGS "${AMR_WIND_CXX_FLAGS}" PARENT_SCOPE)
  set(AMR_WIND_C_FLAGS "${AMR_WIND_C_FLAGS}" PARENT_SCOPE)
  set(AMR_WIND_Fortran_FLAGS "${AMR_WIND_Fortran_FLAGS}" PARENT_SCOPE)

endfunction(set_compile_flags)
