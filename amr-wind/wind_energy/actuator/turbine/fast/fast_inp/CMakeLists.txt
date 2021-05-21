# set variables to make it easy to link/copy this directory
set(NREL5MW_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(NREL5MW_INSTALL_DIR ${CMAKE_BINARY_DIR})
message("NREL Install dir: " ${NREL5MW_INSTALL_DIR})
message("NREL Source dir: " ${NREL5MW_SOURCE_DIR})

file(COPY ${NREL5MW_SOURCE_DIR} DESTINATION ${CMAKE_BINARY_DIR}/)
set(NREL5MW_DISCON_LIB "${CMAKE_SHARED_LIBRARY_PREFIX}DISCON${CMAKE_SHARED_LIBRARY_SUFFIX}")
configure_file(${NREL5MW_INSTALL_DIR}/fast_inp/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ServoDyn.dat.in
               ${NREL5MW_INSTALL_DIR}/fast_inp/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ServoDyn.dat @ONLY)
configure_file(${NREL5MW_INSTALL_DIR}/fast_inp/nrel5mw.fst.in
               ${NREL5MW_INSTALL_DIR}/fast_inp/nrel5mw.fst @ONLY 
               ) 

# Create copies of NREL5MW data and fill correct paths for tests
# This allows turbine files to be unique to each test as opposed to overwriting them
# when we use symbolic links
function(nrel5mw_copy_to_dir copy_to_dir)
    # Set name for library using the prefix and suffix for the current system
    file(COPY ${CMAKE_BINARY_DIR}/fast_inp DESTINATION ${copy_to_dir}/)
    set(NREL5MW_DISCON_LIB "${CMAKE_SHARED_LIBRARY_PREFIX}DISCON${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(NREL5MW_INSTALL_DIR ${copy_to_dir}/fast_inp)
    configure_file(${CMAKE_BINARY_DIR}/fast_inp/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ServoDyn.dat.in
                   ${NREL5MW_INSTALL_DIR}/5MW_Baseline/NRELOffshrBsline5MW_Onshore_ServoDyn.dat @ONLY)
    configure_file(${CMAKE_BINARY_DIR}/fast_inp/nrel5mw.fst.in
                   ${NREL5MW_INSTALL_DIR}/nrel5mw.fst @ONLY 
                   ) 
endfunction(nrel5mw_copy_to_dir) 

add_subdirectory(5MW_Baseline/ServoData)