# Set location of gold files according to system/compiler/compiler_version
set(FCOMPARE_GOLD_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/AMR-WindGoldFiles/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})

if(AMR_WIND_TEST_WITH_FCOMPARE)
  message(STATUS "Test golds directory for fcompare: ${FCOMPARE_GOLD_FILES_DIRECTORY}")
endif()

if(AMR_WIND_ENABLE_MASA AND NOT AMR_WIND_ENABLE_MPI)
  message(WARNING "Running verification tests without MPI enabled will require long run times")
endif()

# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r TEST_NAME)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${FCOMPARE_GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    # Copy files to test working directory
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests in this category
    set(RUNTIME_OPTIONS "time.max_step=10 amr.plot_file=plt time.plot_interval=10 amrex.throw_exception=1 amrex.signal_handling=0")
    if(AMR_WIND_ENABLE_CUDA)
      set(FCOMPARE_TOLERANCE "-r 1e-10")
    endif()
    # Use fcompare to test diffs in plots against gold files
    if(AMR_WIND_TEST_WITH_FCOMPARE)
      set(FCOMPARE_COMMAND "&& ${FCOMPARE_EXE} ${FCOMPARE_TOLERANCE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
    if(AMR_WIND_ENABLE_MPI)
      set(NP 4)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(NP 1)
      unset(MPI_COMMANDS)
    endif()
    # Add test and actual test commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 5400
                         PROCESSORS ${NP}
                         WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                         LABELS "regression"
                         ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_r)

# Regression tests excluded from CI
function(add_test_re TEST_NAME)
    add_test_r(${TEST_NAME})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no_ci")
endfunction(add_test_re)

# Regression test and excluded from CI with dependency
function(add_test_red TEST_NAME TEST_DEPENDENCY)
    add_test_re(${TEST_NAME})
    set_tests_properties(${TEST_NAME} PROPERTIES FIXTURES_REQUIRED fixture_${TEST_DEPENDENCY})
    set_tests_properties(${TEST_DEPENDENCY} PROPERTIES FIXTURES_SETUP fixture_${TEST_DEPENDENCY})
endfunction(add_test_red)

# Verification test using multiple resolutions
function(add_test_v TEST_NAME LIST_OF_GRID_SIZES)
    # Make sure run command is cleared before we construct it
    unset(MASTER_RUN_COMMAND)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Get last item in resolution list so we can find out when we are on the last item in our loop
    list(GET LIST_OF_GRID_SIZES -1 LAST_GRID_SIZE_IN_LIST)
    # Create the commands to run for each resolution
    foreach(GRID_SIZE IN LISTS LIST_OF_GRID_SIZES)
      # Make working directory for test
      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE})
      # Gather all files in source directory for test
      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
      # Copy files to test working directory
      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
      # Set number of cells at runtime according to dimension
      set(NCELLS "${GRID_SIZE} ${GRID_SIZE} ${GRID_SIZE}")
      if(AMR_WIND_ENABLE_MPI)
        set(NP 4)
        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
      else()
        set(NP 1)
        unset(MPI_COMMANDS)
      endif()
      # Set the run command for this resolution
      set(RUN_COMMAND_${GRID_SIZE} "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/${TEST_NAME}.i")
      # Set some runtime options for each resolution
      set(RUNTIME_OPTIONS_${GRID_SIZE} "amrex.throw_exception=1 amrex.signal_handling=0 amr.n_cell=${NCELLS}")
      # Construct our large run command for the entire test with everything &&'d together
      string(APPEND MASTER_RUN_COMMAND "cd ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${RUN_COMMAND_${GRID_SIZE}} ${RUNTIME_OPTIONS_${GRID_SIZE}} > ${TEST_NAME}_${GRID_SIZE}.log")
      # Add another " && " unless we are on the last resolution in the list
      if(NOT ${GRID_SIZE} EQUAL ${LAST_GRID_SIZE_IN_LIST})
        string(APPEND MASTER_RUN_COMMAND " && ")
      endif()
    endforeach()
    # Convert list of grid sizes to space separated string
    list(JOIN LIST_OF_GRID_SIZES " " STRING_OF_GRID_SIZES)
    # Add test and actual test commands to CTest database (python script requires a very specific python environment)
    add_test(${TEST_NAME} sh -c "${MASTER_RUN_COMMAND} && cd ${CURRENT_TEST_BINARY_DIR} && ${PYTHON_EXECUTABLE} ${CURRENT_TEST_SOURCE_DIR}/plotter.py -f ${STRING_OF_GRID_SIZES}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 14400 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}" LABELS "verification;no_ci" ATTACHED_FILES "${CURRENT_TEST_BINARY_DIR}/plots.pdf")
endfunction(add_test_v)

# Standard unit test
function(add_test_u TEST_NAME)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Use a single process
    set(NP 1)
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
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 500
                         PROCESSORS ${NP}
                         WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                         LABELS "unit")
endfunction(add_test_u)

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit_tests)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r(abl_godunov)
add_test_r(boussinesq_bubble_godunov)
add_test_r(tgv_godunov)
add_test_r(abl_mol)

#=============================================================================
# Regression tests excluded from CI
#=============================================================================
add_test_re(tgv_mol)
add_test_re(boussinesq_bubble_mol)
add_test_re(rayleigh_taylor_godunov)
add_test_re(rayleigh_taylor_mol)
add_test_re(abl_godunov_plm)
add_test_re(abl_godunov_explicit)
add_test_re(abl_godunov_cn)
add_test_re(abl_mol_explicit)
add_test_re(abl_mol_cn)
add_test_re(tgv_godunov_plm)
add_test_re(abl_godunov_static_refinement)
add_test_re(abl_ksgsm84_godunov)
add_test_re(abl_godunov_nolim)
add_test_re(ekman_spiral)
add_test_re(vortex_patch_godunov)
add_test_re(zalesak_disk_godunov)

if (NOT AMR_WIND_ENABLE_CUDA)
  add_test_re(ctv_godunov_plm)
endif()

if (AMR_WIND_ENABLE_NETCDF)
  add_test_re(abl_bndry_output)
endif()

if(AMR_WIND_ENABLE_MASA)
  add_test_re(mms_godunov)
  add_test_re(mms_godunov_plm)
  add_test_re(mms_mol)
endif()

#=============================================================================
# Regression tests excluded from CI with a test dependency
#=============================================================================
add_test_red(abl_godunov_restart abl_godunov)
if (AMR_WIND_ENABLE_NETCDF)
  add_test_red(abl_bndry_input abl_bndry_output)
endif()

#=============================================================================
# Verification tests
#=============================================================================
if(AMR_WIND_ENABLE_MASA)
  set(LIST_OF_GRID_SIZES 8 16 32 64)
  add_test_v(mms "${LIST_OF_GRID_SIZES}")
endif()

#=============================================================================
# Performance tests
#=============================================================================

