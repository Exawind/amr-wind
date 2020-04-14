# Set location of gold files according to system/compiler/compiler_version
set(FCOMPARE_GOLD_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AMR-WindGoldFiles/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})

if(AMR_WIND_TEST_WITH_FCOMPARE)
  message(STATUS "Test golds directory for fcompare: ${FCOMPARE_GOLD_FILES_DIRECTORY}")
endif()

# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

# Set TOLERANCE for testing
#if(NOT ${TEST_TOLERANCE} STREQUAL "")
#  set(TOLERANCE ${TEST_TOLERANCE}) # User defined
#else(NOT ${TEST_TOLERANCE} STREQUAL "")
#  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
#        OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
#      set(TOLERANCE "1e-3")
#    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#      set(TOLERANCE "1e-3")
#    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#      set(TOLERANCE "1e-2")
#    else()
#      set(TOLERANCE "1e-8") # Mac default
#    endif()
#  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#      set(TOLERANCE "1e-5")
#    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#      set(TOLERANCE "1e-15")
#    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#      set(TOLERANCE "1e-2")
#    else()
#      set(TOLERANCE "1e-8") # Linux default
#    endif()
#  endif()
#endif()
#message(STATUS "Using test tolerance of ${TOLERANCE}")

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r TEST_NAME NP)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${FCOMPARE_GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Find fcompare
    if(AMR_WIND_TEST_WITH_FCOMPARE)
      set(FCOMPARE ${CMAKE_BINARY_DIR}/${AMREX_SUBMOD_LOCATION}/Tools/Plotfile/fcompare)
    endif()
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    # Copy files to test working directory
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests in this category
    set(RUNTIME_OPTIONS "time.max_step=10 amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1")
    # Use fcompare to test diffs in plots against gold files
    if(AMR_WIND_TEST_WITH_FCOMPARE)
      set(FCOMPARE_COMMAND "&& ${FCOMPARE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
    if(AMR_WIND_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    # Add test and actual test commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 5400 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression")
endfunction(add_test_r)

# Regression tests excluded from CI
function(add_test_re TEST_NAME NP)
    add_test_r(${TEST_NAME} ${NP})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no_ci")
endfunction(add_test_re)

# Standard unit test
function(add_test_u TEST_NAME NP)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    if(AMR_WIND_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    # Add test and commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_unit_test_exe_name}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "unit")
endfunction(add_test_u)

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit_tests 1)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r(tgv_mol 4)
add_test_r(tgv_godunov 4)
add_test_r(boussinesq_bubble_mol 4)
add_test_r(boussinesq_bubble_godunov 4)
add_test_r(abl_mol 4)
add_test_r(abl_godunov 4)

#=============================================================================
# Regression tests excluded from CI
#=============================================================================
add_test_re(rayleigh_taylor_godunov 4)
add_test_re(abl_godunov_explicit 4)
add_test_re(abl_godunov_cn 4)
add_test_re(abl_mol_explicit 4)
add_test_re(abl_mol_cn 4)
add_test_re(tgv_godunov_plm 4)

#=============================================================================
# Verification tests
#=============================================================================

#=============================================================================
# Performance tests
#=============================================================================

