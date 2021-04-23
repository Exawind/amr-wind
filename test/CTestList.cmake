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
macro(setup_test)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(PLOT_GOLD ${FCOMPARE_GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00010)
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    set(RUNTIME_OPTIONS "time.max_step=10 amr.plot_file=plt time.plot_interval=10 amrex.throw_exception=1 amrex.signal_handling=0")
    if(AMR_WIND_ENABLE_MPI)
      if(AMR_WIND_ENABLE_CUDA)
        set(TEST_NP 2)
      else()
        set(TEST_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(TEST_NP 1)
      unset(MPI_COMMANDS)
    endif()
    if(AMR_WIND_ENABLE_CUDA)
      set(FCOMPARE_TOLERANCE "-r 1e-10 --abs_tol 1.0e-12")
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} io.skip_outputs=p")
    endif()
    if(AMR_WIND_TEST_WITH_FCOMPARE)
      set(FCOMPARE_COMMAND "&& ${MPI_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_TOLERANCE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME)
    setup_test()
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 5400
                         PROCESSORS ${TEST_NP}
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
    setup_test()
    unset(MASTER_RUN_COMMAND)
    # Get last item in resolution list so we can find out when we are on the last item in our loop
    list(GET LIST_OF_GRID_SIZES -1 LAST_GRID_SIZE_IN_LIST)
    foreach(GRID_SIZE IN LISTS LIST_OF_GRID_SIZES)
      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE})
      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
      set(NCELLS "${GRID_SIZE} ${GRID_SIZE} ${GRID_SIZE}")
      set(RUN_COMMAND_${GRID_SIZE} "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/${TEST_NAME}.i")
      set(RUNTIME_OPTIONS_${GRID_SIZE} "amrex.throw_exception=1 amrex.signal_handling=0 amr.n_cell=${NCELLS}")
      string(APPEND MASTER_RUN_COMMAND "cd ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${RUN_COMMAND_${GRID_SIZE}} ${RUNTIME_OPTIONS_${GRID_SIZE}} > ${TEST_NAME}_${GRID_SIZE}.log")
      # Add another " && " unless we are on the last resolution in the list
      if(NOT ${GRID_SIZE} EQUAL ${LAST_GRID_SIZE_IN_LIST})
        string(APPEND MASTER_RUN_COMMAND " && ")
      endif()
    endforeach()
    list(JOIN LIST_OF_GRID_SIZES " " STRING_OF_GRID_SIZES)
    add_test(${TEST_NAME} sh -c "${MASTER_RUN_COMMAND} && cd ${CURRENT_TEST_BINARY_DIR} && ${PYTHON_EXECUTABLE} ${CURRENT_TEST_SOURCE_DIR}/plotter.py -f ${STRING_OF_GRID_SIZES}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 14400 PROCESSORS ${TEST_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}" LABELS "verification;no_ci" ATTACHED_FILES "${CURRENT_TEST_BINARY_DIR}/plots.pdf")
endfunction(add_test_v)

# Standard unit test
function(add_test_u TEST_NAME)
    setup_test()
    set(TEST_NP 1)
    if(AMR_WIND_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_unit_test_exe_name}")
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 500
                         PROCESSORS ${TEST_NP}
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
add_test_r(abl_mol)
add_test_r(boussinesq_bubble_godunov)
add_test_r(freestream_godunov)
add_test_r(tgv_godunov)

#=============================================================================
# Regression tests excluded from CI
#=============================================================================
add_test_re(abl_godunov_cn)
add_test_re(abl_godunov_explicit)
add_test_re(abl_godunov_nolim)
add_test_re(abl_godunov_plm)
add_test_re(abl_godunov_static_refinement)
add_test_re(abl_godunov_scalar_velocity_solve)
add_test_re(abl_ksgsm84_godunov)
add_test_re(abl_mol_cn)
add_test_re(abl_mol_explicit)
add_test_re(abl_stable)
add_test_re(abl_unstable)
add_test_re(abl_unstable_constant_wall_model)
add_test_re(abl_unstable_local_wall_model)
add_test_re(abl_unstable_schumann_wall_model)
add_test_re(act_fixed_wing)
add_test_re(act_flat_plate)
add_test_re(boussinesq_bubble_mol)
add_test_re(channel_kwsst)
add_test_re(channel_kwsstiddes)
add_test_re(ekman_spiral)
add_test_re(rayleigh_taylor_godunov)
add_test_re(rayleigh_taylor_mol)
add_test_re(tgv_godunov_plm)
add_test_re(tgv_mol)
add_test_re(vortex_patch_godunov)
add_test_re(zalesak_disk_godunov)
add_test_re(dam_break_godunov)
add_test_re(sloshing_tank)
add_test_re(abl_godunov_weno)

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

# TODO: Enable hypre capability on GPUs
if (AMR_WIND_ENABLE_HYPRE)
  add_test_re(abl_godunov_hypre)
  add_test_re(channel_kwsst_hypre)
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

